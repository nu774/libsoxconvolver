#include <math.h>
#define sin(x) sinf(x)
#define cos(x) cosf(x)
#undef M_PI
/*
 * 0x40490fda = 3.14159250
 * 0x40490fdb = 3.14159274
 */
#define M_PI 3.14159265f
#include "pffft.c"

void pffft_zconvolve(PFFFT_Setup *s, const float *a, const float *b, float *ab) {
  int i, Ncvec = s->Ncvec;
  const v4sf * va = (const v4sf*)a;
  const v4sf * RESTRICT vb = (const v4sf*)b;
  v4sf * vab = (v4sf*)ab;

  float ar, ai, br, bi;

#ifdef __arm__
#error
#endif
  assert(VALIGNED(a) && VALIGNED(b) && VALIGNED(ab));
  ar = ((v4sf_union*)va)[0].f[0];
  ai = ((v4sf_union*)va)[1].f[0];
  br = ((v4sf_union*)vb)[0].f[0];
  bi = ((v4sf_union*)vb)[1].f[0];

  for (i=0; i < Ncvec; i += 2) {
    v4sf ar, ai, br, bi;
    ar = va[2*i+0]; ai = va[2*i+1];
    br = vb[2*i+0]; bi = vb[2*i+1];
    VCPLXMUL(ar, ai, br, bi);
    vab[2*i+0] = ar;
    vab[2*i+1] = ai;
    ar = va[2*i+2]; ai = va[2*i+3];
    br = vb[2*i+2]; bi = vb[2*i+3];
    VCPLXMUL(ar, ai, br, bi);
    vab[2*i+2] = ar;
    vab[2*i+3] = ai;
  }
  if (s->transform == PFFFT_REAL) {
    ((v4sf_union*)vab)[0].f[0] = ar*br;
    ((v4sf_union*)vab)[1].f[0] = ai*bi;
  }
}
