#include <stdio.h>
#include <string.h>
#include "msppriv.h"
#include "ketopt.h"

#define MSP_VERSION "0.0-dirty"

int main_build(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: minisplice <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  bed2bed      convert BED12 to BED\n");
	fprintf(fp, "  version      print the version number\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	int ret = 0;
	if (argc == 1) return usage(stdout);
	else if (strcmp(argv[1], "bed2bed") == 0) return 1;
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
