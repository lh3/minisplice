#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "msppriv.h"
#include "ketopt.h"

#define MSP_VERSION "0.0-dirty"

int main_bed2bed(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: minisplice <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  bed2bed      convert BED12 to BED6\n");
	fprintf(fp, "  version      print the version number\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	int ret = 0;
	msp_realtime();
	if (argc == 1) return usage(stdout);
	else if (strcmp(argv[1], "bed2bed") == 0) main_bed2bed(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		printf("%s\n", MSP_VERSION);
		return 0;
	} else {
		fprintf(stderr, "ERROR: unknown command '%s'\n", argv[1]);
		return 1;
	}

	if (msp_verbose >= 3 && argc > 2 && ret == 0) {
		int i;
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MSP_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, msp_realtime(), msp_cputime(), msp_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}

int main_bed2bed(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	kstring_t out = { 0, 0, 0 };
	msp_bed_t *bed;
	msp_bedview_t *bv = 0;
	int64_t i;
	int c, to_sort = 0, for_train = 0;

	while ((c = ketopt(&o, argc, argv, 1, "st", 0)) >= 0) {
		if (c == 's') to_sort = 1;
		else if (c == 't') for_train = to_sort = 1;
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: minisplice bed2bed [options] <in.bed>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -s     sort\n");
		fprintf(stderr, "  -t     choose training genes\n");
		return 1;
	}
	bed = msp_bed_read(argv[o.ind]);
	if (to_sort) msp_bed_sort(bed);
	if (for_train) bv = msp_bed_select_train(bed);
	if (bv) {
		for (i = 0; i < bv->n; ++i) {
			msp_bed_format(&out, bv->a[i]);
			puts(out.s);
		}
		free(bv->a);
		free(bv);
	} else {
		for (i = 0; i < bed->n; ++i) {
			msp_bed_format(&out, bed->a[i]);
			puts(out.s);
		}
	}
	msp_bed_destroy(bed);
	return 0;
}
