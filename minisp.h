#ifndef MINISP_H
#define MINISP_H

#include <stdint.h>

typedef enum { MSP_FT_FASTX, MSP_FT_BED } msp_ft_t;

typedef struct {
	msp_ft_t type;
	void *fp;
	void *buf;
} msp_file_t;

typedef struct {
	int64_t st, en;
} msp_blk1_t;

typedef struct {
	int64_t st, en;
	int64_t st2, en2; // thickStart and thickEnd
	int64_t score;
	int32_t n_blk, strand;
	char *ctg, *name;
	msp_blk1_t blk[];
} msp_bed1_t;

#endif
