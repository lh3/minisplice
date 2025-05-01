#ifndef MSPPRIV_H
#define MSPPRIV_H

#include "minisp.h"

#define MSP_MALLOC(type, cnt)       ((type*)malloc((cnt) * sizeof(type)))
#define MSP_CALLOC(type, cnt)       ((type*)calloc((cnt), sizeof(type)))
#define MSP_REALLOC(type, ptr, cnt) ((type*)realloc((ptr), (cnt) * sizeof(type)))

#endif
