/* Abstract effect: dft filter
 * Copyright (c) 2008 robs@users.sourceforge.net
 * Copyright (c) 2013 nu774
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#if HAVE_CONFIG_H
#include "config.h"
#endif
#include "soxint.h"
#include "fifo.h"
#include "pffft_ex.h"
#include "fft4g_single.h"

typedef float sample_t;

typedef struct dft_filter_t dft_filter_t;
typedef struct dft_filter_priv_t dft_filter_priv_t;

struct dft_filter_t {
  unsigned dft_length, num_taps, post_peak;
  sample_t * coefs;
  void (* setup)(dft_filter_priv_t *);
  void (* destroy)(dft_filter_t *);
  void (* convolve)(dft_filter_priv_t *, sample_t *, sample_t *);
  void *dft_setup;
};

struct dft_filter_priv_t {
  sample_t       * scratch, * work;
  uint64_t       samples_in, samples_out;
  fifo_t         input_fifo, output_fifo;
  dft_filter_t   filter, * filter_ptr;
};

typedef dft_filter_t filter_t;

static float d2f(double v);

static
void lsx_set_dft_filter(dft_filter_t *f, double *h, int n, int post_peak,
                        int multiplier)
{
  unsigned i;
  f->num_taps = n;
  f->post_peak = post_peak;
  f->dft_length = _lsx_set_dft_length(f->num_taps);
  f->coefs = lsx_aligned_malloc(64, f->dft_length * sizeof(sample_t));
  memset(f->coefs, 0, f->dft_length * sizeof(sample_t));
  for (i = 0; i < f->num_taps; ++i)
    f->coefs[(i + f->dft_length - f->num_taps + 1) & (f->dft_length - 1)]
      = d2f(h[i] / f->dft_length * multiplier);
}

static int start(dft_filter_priv_t * p)
{
  filter_t const * f = p->filter_ptr;
  void * x;
  fifo_create(&p->input_fifo, sizeof(sample_t));
  x = fifo_reserve(&p->input_fifo, f->post_peak);
  memset(x, 0, sizeof(sample_t) * f->post_peak);
  fifo_create(&p->output_fifo, sizeof(sample_t));
  p->scratch = lsx_aligned_malloc(64, f->dft_length * sizeof(sample_t));
  p->work = lsx_aligned_malloc(64, f->dft_length * sizeof(sample_t));
  return 0;
}

static void filter(dft_filter_priv_t * p)
{
  size_t num_in = fifo_occupancy(&p->input_fifo);
  filter_t const * f = p->filter_ptr;
  int const overlap = f->num_taps - 1;
  sample_t * output;

  while (num_in >= f->dft_length) {
    sample_t * input = (sample_t *)fifo_read_ptr(&p->input_fifo);
    fifo_read(&p->input_fifo, f->dft_length - overlap, NULL);
    num_in -= f->dft_length - overlap;
    output = (sample_t *)fifo_reserve(&p->output_fifo, f->dft_length);
    fifo_trim_by(&p->output_fifo, overlap);
    f->convolve(p, input, output);
  }
}

static int flow(dft_filter_priv_t * p, const sample_t * ibuf,
                sample_t * obuf, size_t istride, size_t ostride,
                size_t * isamp, size_t * osamp)
{
  size_t i, odone = min(*osamp, fifo_occupancy(&p->output_fifo));
  sample_t const * s = (sample_t *)fifo_read(&p->output_fifo, odone, NULL);

  for (i = 0; i < odone; ++i, obuf += ostride)
    *obuf = *s++;
  p->samples_out += odone;

  if (*isamp && odone < *osamp) {
    sample_t * t = (sample_t *)fifo_write(&p->input_fifo, *isamp, NULL);
    p->samples_in += *isamp;

    for (i = *isamp; i; --i, ibuf += istride)
      *t++ = *ibuf;
    filter(p);
  }
  else *isamp = 0;
  *osamp = odone;
  return 0;
}

static int drain(dft_filter_priv_t * p, sample_t * obuf, size_t ostride,
                 size_t * osamp)
{
  static size_t isamp = 0;
  uint64_t samples_out = p->samples_in;
  int64_t remaining = samples_out - p->samples_out;
  sample_t * buff = (sample_t *)calloc(1024, sizeof(*buff));

  if (remaining > 0) {
    while ((int64_t)fifo_occupancy(&p->output_fifo) < remaining) {
      fifo_write(&p->input_fifo, 1024, buff);
      p->samples_in += 1024;
      filter(p);
    }
    fifo_trim_to(&p->output_fifo, remaining);
    p->samples_in = 0;
  }
  free(buff);
  return flow(p, 0, obuf, 0, ostride, &isamp, osamp);
}

static int stop(dft_filter_priv_t * p)
{
  filter_t *f = p->filter_ptr;
  fifo_delete(&p->input_fifo);
  fifo_delete(&p->output_fifo);
  lsx_aligned_free(p->scratch);
  lsx_aligned_free(p->work);
  if (f->coefs) {
    lsx_aligned_free(f->coefs);
    f->destroy(f);
    memset(f, 0, sizeof(filter_t));
  }
  return 0;
}



static float d2f(double v)
{
  const float anti_denormal = 1.0e-30f;
  float x = (float)v;
  x += anti_denormal;
  x -= anti_denormal;
  return x;
}

static void dft_setup_pffft(dft_filter_priv_t * p)
{
  filter_t * f = p->filter_ptr;
  f->dft_setup = pffft_new_setup(f->dft_length, PFFFT_REAL);
  pffft_transform(f->dft_setup, f->coefs, f->coefs, p->scratch, PFFFT_FORWARD);
}

static void dft_setup_fft4g(dft_filter_priv_t * p)
{
  filter_t * f = p->filter_ptr;
  f->dft_setup = fft4g_new_setup(f->dft_length);
  fft4g_rdft_forward(f->dft_setup, f->coefs);
}

static void dft_destroy_pffft(dft_filter_t * f)
{
  if (f->dft_setup) pffft_destroy_setup(f->dft_setup);
}

static void dft_destroy_fft4g(dft_filter_t * f)
{
  if (f->dft_setup) fft4g_destroy_setup(f->dft_setup);
}

static void pffft_convolve(dft_filter_priv_t * p, sample_t * s, sample_t * d)
{
  filter_t * f = p->filter_ptr;
  sample_t * w = p->work;

  memcpy(w, s, f->dft_length * sizeof(sample_t));
  pffft_transform(f->dft_setup, w, w, p->scratch, PFFFT_FORWARD);
  pffft_zconvolve(f->dft_setup, w, f->coefs, w);
  pffft_transform(f->dft_setup, w, w, p->scratch, PFFFT_BACKWARD);
  memcpy(d, w, f->dft_length * sizeof(sample_t));
}

static void fft4g_convolve(dft_filter_priv_t * p, sample_t * s, sample_t * d)
{
  unsigned i;
  filter_t * f = p->filter_ptr;

  memcpy(d, s, f->dft_length * sizeof(sample_t));
  fft4g_rdft_forward(f->dft_setup, d);
  d[0] *= f->coefs[0];
  d[1] *= f->coefs[1];
  for (i = 2; i < f->dft_length; i += 2) {
    sample_t tmp = d[i];
    d[i  ] = f->coefs[i  ] * tmp - f->coefs[i+1] * d[i+1];
    d[i+1] = f->coefs[i+1] * tmp + f->coefs[i  ] * d[i+1];
  }
  fft4g_rdft_backword(f->dft_setup, d);
}

typedef struct lsx_convolver_t {
  unsigned nchannels;
  dft_filter_priv_t dft[1];
} lsx_convolver_t;


#if defined(i386) && (defined(__GNUC__) || defined(__clang__))
#include <cpuid.h>
unsigned cpu_feature()
{
  unsigned eax, ebx, ecx, edx;
  __get_cpuid(1, &eax, &ebx, &ecx, &edx);
  return edx;
}
#elif defined(_M_IX86) && defined(_MSC_VER)
#include <intrin.h>
unsigned cpu_feature()
{
  int info[4] = { 0 };
  __cpuid(info, 1);
  return info[3];
}
#endif

lsx_convolver_t *lsx_convolver_create(unsigned nchannels,
                                      double *coefs, unsigned ncoefs,
                                      unsigned post_peak)
{
  unsigned i;
  lsx_convolver_t *state;
  dft_filter_priv_t *p;
  filter_t *f;
  int has_sse = 0;

  if (!nchannels || !coefs || !ncoefs || ncoefs >= 0x10000)
    return 0;
  if ((state = calloc(1, offsetof(lsx_convolver_t, dft[nchannels]))) == 0)
    return 0;
  state->nchannels = nchannels;
  p = &state->dft[0];
  f = &p->filter;
#ifndef DISABLE_SIMD
#if defined(__x86_64__) || defined(_M_X64)
  has_sse = 1;
#elif defined(i386) || defined(_M_IX86)
  has_sse = (cpu_feature() & (1<<25));
#endif
#endif
  lsx_set_dft_filter(f, coefs, ncoefs, post_peak, has_sse ? 1 : 2);
  if (has_sse) {
    f->setup = dft_setup_pffft;
    f->destroy = dft_destroy_pffft;
    f->convolve = pffft_convolve;
  } else {
    f->setup = dft_setup_fft4g;
    f->destroy = dft_destroy_fft4g;
    f->convolve = fft4g_convolve;
  }
  for (i = 0; i < nchannels; ++i) {
    state->dft[i].filter_ptr = f;
    start(&state->dft[i]);
  }
  f->setup(p);
  return state;
}

void lsx_convolver_close(lsx_convolver_t *state)
{
  unsigned n;
  for (n = 0; n < state->nchannels; ++n)
    stop(&state->dft[n]);
  free(state);
}

void lsx_convolver_process(lsx_convolver_t *state, const sample_t *ibuf,
                           sample_t *obuf, size_t *ilen, size_t *olen)
{
  size_t in = *ilen;
  size_t out = *olen;
  unsigned i, nc = state->nchannels;

  for (i = 0; i < nc; ++i) {
    *ilen = in;
    *olen = out;
    if (in > 0)
      flow(&state->dft[i], ibuf + i, obuf + i, nc, nc, ilen, olen);
    else
      drain(&state->dft[i], obuf + i, nc, olen);
  }
}

void lsx_convolver_process_ni(lsx_convolver_t *state,
                              const sample_t * const *ibuf, sample_t **obuf,
                              size_t istride, size_t ostride,
                              size_t *ilen, size_t *olen)
{
  size_t in = *ilen;
  size_t out = *olen;
  unsigned i, nc = state->nchannels;

  for (i = 0; i < nc; ++i) {
    *ilen = in;
    *olen = out;
    if (in > 0)
      flow(&state->dft[i], ibuf[i], obuf[i], istride, ostride, ilen, olen);
    else
      drain(&state->dft[i], obuf[i], ostride, olen);
  }
}

#include "version.h"
const char *lsx_convolver_version_string(void)
{
  return SOXCONVOLVER_VERSION_STRING;
}

void lsx_free(void *memory)
{
  free(memory);
}
