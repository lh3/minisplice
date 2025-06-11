#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "msppriv.h"
#include "ketopt.h"

int main_bed2bed(int argc, char *argv[]);
int main_gentrain(int argc, char *argv[]);
int main_train(int argc, char *argv[]);
int main_inspect(int argc, char *argv[]);
int main_predict(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: minisplice <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  gentrain     generate training data\n");
	fprintf(fp, "  train        model training\n");
	fprintf(fp, "  inspect      print the model structure\n");
	fprintf(fp, "  predict      predict splice sites\n");
//	fprintf(fp, "  bed2bed      convert BED12 to BED6\n");
	fprintf(fp, "  version      print the version number\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	int ret = 0;
	msp_realtime();
	if (argc == 1) return usage(stdout);
	else if (strcmp(argv[1], "bed2bed") == 0) ret = main_bed2bed(argc-1, argv+1);
	else if (strcmp(argv[1], "gentrain") == 0) ret = main_gentrain(argc-1, argv+1);
	else if (strcmp(argv[1], "train") == 0) ret = main_train(argc-1, argv+1);
	else if (strcmp(argv[1], "inspect") == 0) ret = main_inspect(argc-1, argv+1);
	else if (strcmp(argv[1], "predict") == 0) ret = main_predict(argc-1, argv+1);
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

int main_gentrain(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	msp_bed_t *bed;
	msp_file_t *fx;
	int32_t c, ext = 100;
	double frac_pos = 0.25;
	msp_tdata_t *d;

	while ((c = ketopt(&o, argc, argv, 1, "l:p:", 0)) >= 0) {
		if (c == 'l') ext = atoi(o.arg);
		else if (c == 'p') frac_pos = atof(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: minisplice gentrain [options] <in.bed> <in.fastx>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l INT       length of flanking sequences [%d]\n", ext);
		fprintf(stderr, "  -p FLOAT     fraction of positive sites [%g]\n", frac_pos);
		return 1;
	}
	bed = msp_bed_read(argv[o.ind]);
	if (bed == 0 && msp_verbose >= 1) {
		fprintf(stderr, "ERROR: fail to open BED file '%s'\n", argv[o.ind]);
		return 1;
	}
	msp_bed_idxctg(bed);
	fx = msp_fastx_open(argv[o.ind+1]);
	if (bed == 0 && msp_verbose >= 1) {
		msp_bed_destroy(bed);
		fprintf(stderr, "ERROR: fail to open FASTX file '%s'\n", argv[o.ind+1]);
		return 1;
	}
	d = msp_gen_train(bed, fx, ext, frac_pos);
	msp_file_close(fx);
	msp_tdata_dump(stdout, bed, d);
	msp_tdata_destroy(d);
	msp_bed_destroy(bed);
	return 0;
}

int main_bed2bed(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	kstring_t out = { 0, 0, 0 };
	msp_bed_t *bed;
	int64_t i;
	int c, to_sort = 0, out_neg = 0;

	while ((c = ketopt(&o, argc, argv, 1, "sn", 0)) >= 0) {
		if (c == 's') to_sort = 1;
		else if (c == 'n') out_neg = 1;
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: minisplice bed2bed [options] <in.bed>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -n     generate negative regions\n");
		fprintf(stderr, "  -s     sort\n");
		return 1;
	}
	bed = msp_bed_read(argv[o.ind]);
	if (to_sort) msp_bed_sort(bed);
	if (out_neg) {
		int32_t cid;
		msp_bed_idxctg(bed);
		for (cid = 0; cid < bed->h->n; ++cid) {
			int64_t n;
			msp192_t *a;
			a = msp_bed_gen_negreg(bed, cid, &n);
			for (i = 0; i < n; ++i) {
				out.l = 0;
				msp_sprintf_lite(&out, "%s\t%ld\t%ld\t.\t.\t%c\n", bed->h->a[cid], (long)a[i].x, (long)a[i].y, "+-"[a[i].z]);
				fwrite(out.s, 1, out.l, stdout);
			}
			free(a);
		}
	} else {
		for (i = 0; i < bed->n; ++i) {
			msp_bed_format(&out, bed->a[i]);
			puts(out.s);
		}
	}
	msp_bed_destroy(bed);
	free(out.s);
	return 0;
}
