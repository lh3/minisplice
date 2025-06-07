#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "kann.h"
#include "msppriv.h"
#include "ketopt.h"

static kann_t *msp_model_gen(int n_out, int len, int k_size, int n_flt, int n_fc, float dropout)
{
	kad_node_t *t;
	t = kad_feed(3, 1, 4, len), t->ext_flag |= KANN_F_IN;
	t = kad_relu(kann_layer_conv1d(t, n_flt, k_size, 1, 0));
	t = kad_max1d(t, 3, 3, 0);
	t = kad_relu(kann_layer_conv1d(t, n_flt, k_size, 1, 0));
	if (dropout > 0.0f) t = kann_layer_dropout(t, dropout);
	t = kad_max1d(t, 3, 3, 0);
	t = kad_relu(kann_layer_dense(t, n_fc));
	if (dropout > 0.0f) t = kann_layer_dropout(t, dropout);
	return kann_new(kann_layer_cost(t, n_out, KANN_C_CEB), 0);
}

typedef struct {
	int64_t n;
	int32_t n_label, len;
	float **x, **y;
} msp_fdata_t;

void msp_seq2vec(int32_t len, const uint8_t *s, float *x)
{
	int32_t i, c;
	memset(x, 0, len * 4 * sizeof(float));
	for (c = 0; c < 4; ++c) {
		float *x1 = &x[c * len];
		for (i = 0; i < len; ++i)
			if (s[i] == c)
				x1[i] = 1.0f;
	}
}

static msp_fdata_t *msp_s2fdata(const msp_sdata_t *d)
{
	int64_t i;
	msp_fdata_t *f;
	f = MSP_CALLOC(msp_fdata_t, 1);
	f->n = d->n, f->n_label = d->n_label, f->len = d->len;
	f->x = MSP_CALLOC(float*, f->n);
	f->y = MSP_CALLOC(float*, f->n);
	for (i = 0; i < f->n; ++i) {
		f->x[i] = MSP_CALLOC(float, 4 * f->len);
		f->y[i] = MSP_CALLOC(float, f->n_label);
		assert(d->a[i].label < f->n_label);
		f->y[i][d->a[i].label] = 1.0f;
		msp_seq2vec(f->len, d->a[i].seq, f->x[i]);
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

int main_train(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, k_size = 5, n_flt = 16, n_fc = 16, min_epoch = 3, max_epoch = 100, mb_sz = 64, n_thread = 1;
	int max_drop_streak = 10, seed = 11, print_model = 0;
	float lr = 0.001f, dropout = 0.0f;
	msp_sdata_t *d;
	msp_fdata_t *f;
	char *fn_in = 0, *fn_out = 0;
	kann_t *ann;

	while ((c = ketopt(&o, argc, argv, 1, "t:k:f:m:e:E:r:d:s:i:o:pF:", 0)) >= 0) {
		if (c == 't') n_thread = atoi(o.arg);
		else if (c == 'k') k_size = atoi(o.arg);
		else if (c == 'f') n_flt = atoi(o.arg);
		else if (c == 'F') n_fc = atoi(o.arg);
		else if (c == 'E') max_epoch = atoi(o.arg);
		else if (c == 'm') mb_sz = atoi(o.arg);
		else if (c == 'r') lr = atof(o.arg);
		else if (c == 'd') dropout = atof(o.arg);
		else if (c == 's') seed = atoi(o.arg);
		else if (c == 'i') fn_in = o.arg;
		else if (c == 'o') fn_out = o.arg;
		else if (c == 'p') print_model = 1;
		else if (c == 'e') min_epoch = atoi(o.arg);
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: minisplice train [options] <in.data>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  Model construction:\n");
		fprintf(stderr, "    -k INT     1D-CNN kernel size [%d]\n", k_size);
		fprintf(stderr, "    -f INT     number of features per 1D-CNN layer [%d]\n", n_flt);
		fprintf(stderr, "    -F INT     number of neurons in the dense layer [%d]\n", n_fc);
		fprintf(stderr, "    -d FLOAT   dropout rate (use for large models) [%g]\n", dropout);
		fprintf(stderr, "  Model training:\n");
		fprintf(stderr, "    -r FLOAT   learning rate [%g]\n", lr);
		fprintf(stderr, "    -e INT     min number of epoches [%d]\n", min_epoch);
		fprintf(stderr, "    -E INT     max number of epoches [%d]\n", max_epoch);
		fprintf(stderr, "    -m INT     minibatch size [%d]\n", mb_sz);
		fprintf(stderr, "    -s INT     random seed [%d]\n", seed);
		fprintf(stderr, "    -t INT     number of threads [%d]\n", n_thread);
		fprintf(stderr, "  Model I/O:\n");
		fprintf(stderr, "    -i FILE    input model []\n");
		fprintf(stderr, "    -o FILE    output model []\n");
		fprintf(stderr, "    -p         print model structure\n");
		fprintf(stderr, "Input format:\n");
		fprintf(stderr, "  Two columns: integer label, and fixed-length sequence\n");
		fprintf(stderr, "  Or gentrain output\n");
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

	if (fn_in) ann = kann_load(fn_in);
	ann = msp_model_gen(f->n_label, f->len, k_size, n_flt, n_fc, dropout);
	assert(ann);

	if (print_model) {
		kad_print_graph(stdout, ann->n, ann->v);
		goto end_train;
	}

	kann_srand(seed);
	if (n_thread > 1) kann_mt(ann, n_thread, mb_sz);
	kann_train_fnn1b(ann, lr, mb_sz, max_epoch, min_epoch, max_drop_streak, 0.2f, f->n, f->x, f->y);
	if (fn_out) kann_save(fn_out, ann);

end_train:
	kann_delete(ann);
	msp_fdata_destroy(f);
	return 0;
}
