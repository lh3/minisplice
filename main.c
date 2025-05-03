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
	msp_file_t *fp;
	kstring_t out = { 0, 0, 0 };
	uint32_t err;
	msp_bed1_t *p;

	if (argc == 1) {
		fprintf(stderr, "Usage: minisplice bed2bed <in.bed>\n");
		return 1;
	}

	fp = msp_bed_open(argv[1]);
	while ((p = msp_bed_read1(fp, &err)) != NULL) {
		int32_t i;
		out.l = 0;
		for (i = 0; i < p->n_blk; ++i) {
			msp_sprintf_lite(&out, "%s\t%ld\t%ld\t", p->ctg, p->blk[i].st, p->blk[i].en);
			if (p->name) msp_sprintf_lite(&out, "%s", p->name);
			else msp_sprintf_lite(&out, ".");
			if (p->score >= 0) msp_sprintf_lite(&out, "\t%ld", p->score);
			else msp_sprintf_lite(&out, ".");
			msp_sprintf_lite(&out, "\t%c", p->strand == 0? '.' : p->strand > 0? '+' : '-');
			msp_sprintf_lite(&out, "\n");
		}
		fwrite(out.s, 1, out.l, stdout);
		free(p);
	}
	msp_file_close(fp);
	return 0;
}
