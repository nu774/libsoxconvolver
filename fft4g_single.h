#ifndef FFT4G_SINGLE_H
#define FFT4G_SINGLE_H

#ifdef __cplusplus
extern "C" {
#endif

/*
Fast Fourier/Cosine/Sine Transform
    dimension   :one
    data length :power of 2
    decimation  :frequency
    radix       :4, 2
    data        :inplace
    table       :use
functions
    cdft: Complex Discrete Fourier Transform
    rdft: Real Discrete Fourier Transform
    ddct: Discrete Cosine Transform
    ddst: Discrete Sine Transform
    dfct: Cosine Transform of RDFT (Real Symmetric DFT)
    dfst: Sine Transform of RDFT (Real Anti-symmetric DFT)
*/
void cdft(int n, int isgn, float *a, int *ip, float *w);
void rdft(int n, int isgn, float *a, int *ip, float *w);
void ddct(int n, int isgn, float *a, int *ip, float *w);
void ddst(int n, int isgn, float *a, int *ip, float *w);
void dfct(int n, float *a, float *t, int *ip, float *w);
void dfst(int n, float *a, float *t, int *ip, float *w);

typedef struct fft4g_setup {
  int n;
  int *br;
  float *sc;
} fft4g_setup;

fft4g_setup *fft4g_new_setup(int n);

void fft4g_destroy_setup(fft4g_setup *setup);

static inline void fft4g_rdft_forward(fft4g_setup *p, float *a)
{
    rdft(p->n, 1, a, p->br, p->sc);
}
static inline void fft4g_rdft_backword(fft4g_setup *p, float *a)
{
    rdft(p->n, -1, a, p->br, p->sc);
}

#ifdef __cplusplus
}
#endif

#endif /* FFT4G_SINGLE_H */
