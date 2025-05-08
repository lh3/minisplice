#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "kann.h"
#include "msppriv.h"
#include "ketopt.h"
#include "ksort.h"

#define key_pdata(a) ((a).x)
KRADIX_SORT_INIT(msp_pdata, msp_pdata1_t, key_pdata, 64)

void msp_predict1(msp_pdata_t *t, kann_t *ann, int64_t len, const uint8_t *seq, int32_t mb_sz, int32_t type)
{
	int32_t l, k, ext, alen, n_in, n_out, i_out;
	int64_t i, *mb2i, n0 = t->n;
	uint8_t x, *s;
	float *x1;

	for (i = 0, l = k = 0, x = 0; i < len; ++i) {
		int c = seq[i];
		uint64_t z = (uint64_t)-1;
		if (c < 4) {
			x = (x << 2 | c) & 0xf;
			if (++l >= 2) {
				if (type == 0) { // donor
					if (x == (2<<2|3)) z = (i-1)<<3 | 0<<2 | 0<<1; // GT
					else if (x == (0<<2|1)) z = (i+1)<<3 | 1<<2 | 0<<1; // AC
				} else { // acceptor
					if (x == (0<<2|2)) z = (i+1)<<3 | 0<<2 | 1<<1; // AG
					else if (x == (1<<2|3)) z = (i-1)<<3 | 1<<2 | 1<<1; // CT
				}
				if (z != (uint64_t)-1 && (z>>1&1) == type) {
					MSP_GROW(msp_pdata1_t, t->a, t->n, t->m);
					t->a[t->n++].x = z;
				}
			}
		} else l = 0, x = 0;
	}
	radix_sort_msp_pdata(t->a + n0, t->a + t->n);

	kann_switch(ann, 0);
	i_out = kann_find(ann, KANN_F_OUT, 0);
	n_in = kann_dim_in(ann);
	n_out = kann_dim_out(ann);
	assert(n_in % 4 == 0 && n_out == 2);
	alen = n_in / 4;
	assert(alen % 2 == 0);
	ext = (alen - 2) / 2;

	s = MSP_CALLOC(uint8_t, alen);
	x1 = MSP_CALLOC(float, mb_sz * n_in);
	mb2i = MSP_CALLOC(int64_t, mb_sz);
	kann_feed_bind(ann, KANN_F_IN, 0, &x1);
	for (i = 0; i < t->n; ++i)
		t->a[i].f = -1.0f, t->a[i].s = -1;
	for (i = 0; i < t->n;) {
		int64_t j;
		int32_t k, l;
		const float *y1;
		for (j = i, k = 0; j < t->n && k < mb_sz; ++j) {
			int rc;
			rc = msp_get_seq_in_place(s, len, seq, t->a[j].x, ext);
			if (rc < 0) continue;
			msp_seq2vec(alen, s, &x1[k * n_in]);
			mb2i[k++] = j;
		}
		kann_set_batch_size(ann, k);
		kann_eval_out(ann);
		y1 = ann->v[i_out]->x;
		for (l = 0; l < k; ++l)
			t->a[mb2i[l]].f = y1[l * n_out + 1];
		#if 0 // for debugging only
		for (l = 0; l < k; ++l) {
			int32_t j;
			msp_pdata1_t *p = &t->a[mb2i[l]];
			msp_get_seq_in_place(s, len, seq, p->x, ext);
			msp_seq2vec(alen, s, x1);
			for (j = 0; j < alen; ++j) s[j] = "ACGTN"[s[j]];
			fprintf(stderr, "%ld\t%c\t%c\t%f\t%f\t", (long)(p->x>>3), "+-"[p->x>>2&1], "DA"[p->x>>1&1], y1[l*n_out], y1[l*n_out+1]);
			fwrite(s, 1, alen, stderr);
			fputc('\n', stderr);
		}
		#endif
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

void msp_eval_update(msp_eval_t *e)
{
	int32_t b;
	int64_t tot_n = 0, tot_p = 0;
	e->tot_t = e->tot_p = 0;
	for (b = e->n_bin - 1; b >= 0; --b) {
		e->tot_t += e->bin[b].mt;
		e->tot_p += e->bin[b].mp;
	}
	for (b = e->n_bin - 1; b >= 0; --b) {
		tot_p += e->bin[b].mp;
		tot_n += e->bin[b].mt - e->bin[b].mp;
		e->bin[b].tp = tot_p;
		e->bin[b].fn = e->tot_p - tot_p;
		e->bin[b].fp = tot_n;
		e->bin[b].tn = e->tot_t - e->tot_p - tot_n;
	}
}

msp_eval_t *msp_eval(kann_t *ann, msp_file_t *fx, const msp_bed_t *bed, int32_t mb_sz, int32_t type, float step)
{
	int32_t len, n_bin;
	const char *name, *seq;
	msp_pdata_t u = {0,0,0};
	msp_eval_t *e;
	n_bin = (int32_t)((1.0 + step - 1e-3) / step);
	e = (msp_eval_t*)calloc(1, sizeof(msp_eval_t) + n_bin * sizeof(msp_evalbin_t));
	e->n_bin = n_bin, e->step = step;
	while ((len = msp_fastx_read(fx, &name, &seq)) >= 0) {
		int32_t cid;
		int64_t i, j;
		msp64a_t td = {0,0,0}, ta = {0,0,0}, *pt;

		cid = msp_strmap_get(bed->h, name);
		if (cid < 0) continue;
		msp_gen_pos(&td, &ta, bed, cid, len, (uint8_t*)seq);
		pt = type == 0? &td : &ta;
		u.n = 0;
		msp_predict1(&u, ann, len, (uint8_t*)seq, mb_sz, type);

		i = j = 0;
		while (i < u.n && j < pt->n) {
			if (u.a[i].x>>1 == pt->a[j]>>1) u.a[i].x |= 1, ++i, ++j;
			else if (u.a[i].x>>1 < pt->a[j]>>1) ++i;
			else ++j;
		}
		for (i = 0; i < u.n; ++i) {
			int32_t b = u.a[i].f >= 0.0f? (int32_t)(u.a[i].f / step) : 0;
			if (b >= n_bin) b = n_bin - 1;
			e->bin[b].mt++;
			if (u.a[i].x&1) e->bin[b].mp++;
		}
		free(td.a); free(ta.a);
	}
	free(u.a);
	msp_eval_update(e);
	return e;
}

int main_predict(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t c, n_thread = 1, mb_sz = 128, type = 0; // 0 for donor and 1 for acceptor
	msp_file_t *fx;
	kann_t *ann;
	char *fn_bed = 0;
	float step = 0.05f;

	while ((c = ketopt(&o, argc, argv, 1, "adt:b:e:s:", 0)) >= 0) {
		if (c == 't') n_thread = atoi(o.arg);
		else if (c == 'b') mb_sz = atoi(o.arg);
		else if (c == 'a') type = 1;
		else if (c == 'd') type = 0;
		else if (c == 'e') fn_bed = o.arg;
		else if (c == 's') step = atof(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: minisplice predict [options] <in.kan> <in.fastx>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  Prediction:\n");
		fprintf(stderr, "    -a          acceptor model (donor model by default)\n");
		fprintf(stderr, "    -t INT      number of threads [%d]\n", n_thread);
		fprintf(stderr, "    -b INT      minibatch size [%d]\n", mb_sz);
		fprintf(stderr, "  Evaluation:\n");
		fprintf(stderr, "    -e FILE     annotated splice sites in BED12 []\n");
		fprintf(stderr, "    -s FLOAT    step [%g]\n", step);
		return 1;
	}

	ann = kann_load(argv[o.ind]);
	if (n_thread > 1) kann_mt(ann, n_thread, mb_sz);
	fx = msp_fastx_open(argv[o.ind+1]);
	assert(ann && fx);
	if (fn_bed) {
		int32_t i;
		FILE *fp = stdout;
		msp_bed_t *bed;
		msp_eval_t *e;
		bed = msp_bed_read(fn_bed);
		msp_bed_idxctg(bed);
		assert(bed);
		e = msp_eval(ann, fx, bed, mb_sz, type, step);
		msp_bed_destroy(bed);
		fprintf(fp, "CC\tTT  #allSites\n");
		fprintf(fp, "CC\tTP  #posSites\n");
		fprintf(fp, "CC\tST  step\n");
		fprintf(fp, "CC\tNB  #bins\n");
		fprintf(fp, "CC\tBN  bin  all  pos  TP  FP  TN  FN  FPR  SN\n");
		fprintf(fp, "//\n");
		fprintf(fp, "TT\t%ld\n", (long)e->tot_t);
		fprintf(fp, "TP\t%ld\n", (long)e->tot_p);
		fprintf(fp, "ST\t%g\n", step);
		fprintf(fp, "NB\t%d\n", e->n_bin);
		for (i = e->n_bin - 1; i > 0; --i) {
			const msp_evalbin_t *b = &e->bin[i];
			fprintf(fp, "BN\t%d\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%.4g\t%.4g\n", i, (long)b->mt, (long)b->mp, (long)b->tp,
				(long)b->fp, (long)b->tn, (long)b->fn, (double)b->tp / (b->tp + b->fn), (double)b->fp / (b->fp + b->tn));
		}
		free(e);
	} else {
		msp_predict_print(ann, fx, mb_sz, type);
	}

	msp_file_close(fx);
	kann_delete(ann);
	return 0;
}
