#ifndef PFFFT_EX_H
#define PFFFT_EX_H

#include "pffft.h"

#ifdef __cplusplus
extern "C" {
#endif

void pffft_zconvolve(PFFFT_Setup *s, const float *a, const float *b, float *ab);

#ifdef __cplusplus
}
#endif

#endif /* PFFFT_EX_H */
