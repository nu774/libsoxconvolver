#include <stdlib.h>
#include <math.h>
#include "fft4g_single.h"
#define double float
#define sin(x) sinf(x)
#define cos(x) cosf(x)
#define atan(x) atanf(x)
#include "fft4g.c"
#undef double

#define fft4g_br_len(n) (2 + (1 << (int)(log(n / 2 + .5) / log(2.)) / 2))
#define fft4g_sc_len(n) (n * 5 / 4)

fft4g_setup *fft4g_new_setup(int n)
{
    fft4g_setup *p = malloc(sizeof(fft4g_setup));
    p->n = n;
    p->br = calloc(fft4g_br_len(n), sizeof(int));
    p->sc = calloc(fft4g_sc_len(n), sizeof(float));
    return p;
}

void fft4g_destroy_setup(fft4g_setup *setup)
{
    free(setup->br);
    free(setup->sc);
    free(setup);
}
