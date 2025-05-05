#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "kann.h"
#include "msppriv.h"
#include "ketopt.h"

static kann_t *msp_model_gen(int n_layer, int len, int k_size, int n_flt, int n_fc, float dropout)
{
	kad_node_t *t;
	int i;
	t = kad_feed(3, 1, 4, len), t->ext_flag |= KANN_F_IN;
	for (i = 0; i < n_layer; ++i) {
		t = kad_relu(kann_layer_conv1d(t, n_flt, k_size, 1, 0));
		if (dropout > 0.0f) t = kann_layer_dropout(t, dropout);
		t = kad_max1d(t, 3, 3, 0);
	}
	t = kad_relu(kann_layer_dense(t, n_fc));
	if (dropout > 0.0f) t = kann_layer_dropout(t, dropout);
	t = kad_relu(kann_layer_dense(t, n_fc));
	if (dropout > 0.0f) t = kann_layer_dropout(t, dropout);
	return kann_new(kann_layer_cost(t, 2, KANN_C_CEB), 0);
}

static kann_t *msp_model_gen2(int n_layer, int len, int k_size, int n_flt, int n_fc, float dropout)
{
	kad_node_t *t, *s[3];
	int i, k, len_end = len / 4;
	t = kad_feed(3, 1, 4, len), t->ext_flag |= KANN_F_IN;
	s[0] = kad_slice(t, 2, 0, len_end);
	s[1] = kad_slice(t, 2, len_end, len - len_end);
	s[2] = kad_slice(t, 2, len - len_end, len);
	for (i = 0; i < n_layer; ++i) {
		for (k = 0; k < 3; ++k) {
			s[k] = kad_relu(kann_layer_conv1d(s[k], n_flt, k_size, 1, 0));
			if (dropout > 0.0f) s[k] = kann_layer_dropout(s[k], dropout);
			s[k] = kad_max1d(s[k], 3, 3, 0);
		}
	}
	t = kad_concat(2, 3, s[0], s[1], s[2]);
	t = kad_relu(kann_layer_dense(t, n_fc));
	if (dropout > 0.0f) t = kann_layer_dropout(t, dropout);
	t = kad_relu(kann_layer_dense(t, n_fc));
	if (dropout > 0.0f) t = kann_layer_dropout(t, dropout);
	return kann_new(kann_layer_cost(t, 2, KANN_C_CEB), 0);
}

typedef struct {
	int64_t n;
	int32_t n_label, len;
	float **x, **y;
} msp_fdata_t;

static msp_fdata_t *msp_s2fdata(const msp_sdata_t *d)
{
	int64_t i;
	msp_fdata_t *f;
	f = MSP_CALLOC(msp_fdata_t, 1);
	f->n = d->n, f->n_label = d->n_label, f->len = d->len;
	f->x = MSP_CALLOC(float*, f->n);
	f->y = MSP_CALLOC(float*, f->n);
	for (i = 0; i < f->n; ++i) {
		int32_t j, c;
		f->x[i] = MSP_CALLOC(float, 4 * f->len);
		f->y[i] = MSP_CALLOC(float, f->n_label);
		assert(d->a[i].label < f->n_label);
		f->y[i][d->a[i].label] = 1.0f;
		for (c = 0; c < 4; ++c) {
			float *x1 = &f->x[i][c * f->len];
			for (j = 0; j < f->len; ++j)
				if (d->a[i].seq[j] == c)
					x1[j] = 1.0f;
		}
	}
	return f;
}

static void msp_fdata_destroy(msp_fdata_t *f)
{
	int64_t i;
	for (i = 0; i < f->n; ++i) {
		free(f->x[i]); free(f->y[i]);
	}
	free(f->x); free(f->y);
	free(f);
}

int main_train0(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, n_layer = 2, k_size = 5, n_flt = 32, n_fc = 64, max_epoch = 100, mb_sz = 64, n_thread = 1;
	int max_drop_streak = 10, seed = 11, use_3piece = 0;
	float lr = 0.001f, dropout = 0.2f;
	msp_sdata_t *d;
	msp_fdata_t *f;
	kann_t *ann;

	while ((c = ketopt(&o, argc, argv, 1, "t:k:l:f:m:b:r:d:s:3", 0)) >= 0) {
		if (c == 't') n_thread = atoi(o.arg);
		else if (c == 'k') k_size = atoi(o.arg);
		else if (c == 'l') n_layer = atoi(o.arg);
		else if (c == 'f') n_flt = atoi(o.arg);
		else if (c == 'm') max_epoch = atoi(o.arg);
		else if (c == 'b') mb_sz = atoi(o.arg);
		else if (c == 'r') lr = atof(o.arg);
		else if (c == 'd') dropout = atof(o.arg);
		else if (c == 's') seed = atoi(o.arg);
		else if (c == '3') use_3piece = 1;
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: minisplice train0 [options] <in.data>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l INT     number of 1d CNN layers [%d]\n", n_layer);
		fprintf(stderr, "  -k INT     kernel size [%d]\n", k_size);
		fprintf(stderr, "  -f INT     number of filters [%d]\n", n_flt);
		fprintf(stderr, "  -3         use a 3-piece model\n");
		fprintf(stderr, "  -r FLOAT   learning rate [%g]\n", lr);
		fprintf(stderr, "  -m INT     max number of epoches [%d]\n", max_epoch);
		fprintf(stderr, "  -b INT     minibatch size [%d]\n", mb_sz);
		fprintf(stderr, "  -d FLOAT   dropout [%g]\n", dropout);
		fprintf(stderr, "  -s INT     random seed [%d]\n", seed);
		fprintf(stderr, "  -t INT     number of threads [%d]\n", n_thread);
		return 1;
	}
	d = msp_sdata_read(argv[o.ind]);
	if (d == 0) {
		if (msp_verbose >= 1)
			fprintf(stderr, "ERROR: failed to read data file '%s'\n", argv[o.ind]);
		return 1;
	}
	if (msp_verbose >= 3)
		fprintf(stderr, "[M::%s] read %d labels and %ld sequences of length %d\n", __func__, d->n_label, (long)d->n, d->len);
	f = msp_s2fdata(d);
	msp_sdata_destroy(d);
	kann_srand(seed);
	if (use_3piece) ann = msp_model_gen2(n_layer, f->len, k_size, n_flt, n_fc, dropout);
	else ann = msp_model_gen(n_layer, f->len, k_size, n_flt, n_fc, dropout);
	assert(ann);
	if (n_thread > 1) kann_mt(ann, n_thread, mb_sz);
	kann_train_fnn1(ann, lr, mb_sz, max_epoch, max_drop_streak, 0.2f, f->n, f->x, f->y);
	kann_delete(ann);
	msp_fdata_destroy(f);
	return 0;
}
