#include <stdlib.h>
#include <stdio.h>
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
		t = kad_max1d(t, 3, 3, 0);
	}
	t = kann_layer_dropout(t, dropout);
	t = kann_layer_dense(t, n_fc);
	t = kad_relu(t);
	t = kann_layer_dropout(t, dropout);
	return kann_new(kann_layer_cost(t, 2, KANN_C_CEB), 0);
}

int main_train0(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int n_layer = 2, len = -1, k_size = 5, n_fc = 64, max_epoch = 100, mb_sz = 64, n_thread = 1;
	float dropout = 0.2f;
	int32_t c;
	msp_tdata_t *d;

	while ((c = ketopt(&o, argc, argv, 1, "t:k:", 0)) >= 0) {
		if (c == 't') n_thread = atoi(o.arg);
		else if (c == 'k') k_size = atoi(o.arg);
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: minisplice train0 [options] <in.data>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     kernel size [%d]\n", k_size);
		fprintf(stderr, "  -t INT     number of threads [%d]\n", n_thread);
		return 1;
	}
	d = msp_tdata_read(argv[o.ind]);
	msp_tdata_dump(stdout, 0, d);
	msp_tdata_destroy(d);
	return 0;
}
