/* libSoX internal DSP functions.
 * All public functions & data are prefixed with lsx_ .
 *
 * Copyright (c) 2008 robs@users.sourceforge.net
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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "soxint.h"

static
double lsx_bessel_I_0(double x)
{
  double term = 1, sum = 1, last_sum, x2 = x / 2;
  int i = 1;
  do {
    double y = x2 / i++;
    last_sum = sum, sum += term *= y * y;
  } while (sum != last_sum);
  return sum;
}

int _lsx_set_dft_length(int num_taps) /* Set to 4 x nearest power of 2 */
{      /* or half of that if danger of causing too many cache misses. */
  int min = 10 /* sox_globals.log2_dft_min_size */;
  double d = log((double)num_taps) / log(2.);
  return 1 << range_limit((int)(d + 2.77), min, max((int)(d + 1.77), 17));
}

static
double lsx_kaiser_beta(double att, double tr_bw)
{
  if (att >= 60) {
    static const double coefs[][4] = {
      {-6.784957e-10,1.02856e-05,0.1087556,-0.8988365+.001},
      {-6.897885e-10,1.027433e-05,0.10876,-0.8994658+.002},
      {-1.000683e-09,1.030092e-05,0.1087677,-0.9007898+.003},
      {-3.654474e-10,1.040631e-05,0.1087085,-0.8977766+.006},
      {8.106988e-09,6.983091e-06,0.1091387,-0.9172048+.015},
      {9.519571e-09,7.272678e-06,0.1090068,-0.9140768+.025},
      {-5.626821e-09,1.342186e-05,0.1083999,-0.9065452+.05},
      {-9.965946e-08,5.073548e-05,0.1040967,-0.7672778+.085},
      {1.604808e-07,-5.856462e-05,0.1185998,-1.34824+.1},
      {-1.511964e-07,6.363034e-05,0.1064627,-0.9876665+.18},
    };
    double realm = log(tr_bw/.0005)/log(2.);
    double const * c0 = coefs[range_limit(  (int)realm, 0, (int)array_length(coefs)-1)];
    double const * c1 = coefs[range_limit(1+(int)realm, 0, (int)array_length(coefs)-1)];
    double b0 = ((c0[0]*att + c0[1])*att + c0[2])*att + c0[3];
    double b1 = ((c1[0]*att + c1[1])*att + c1[2])*att + c1[3];
    return b0 + (b1 - b0) * (realm - (int)realm);
  }
  if (att > 50   ) return .1102 * (att - 8.7);
  if (att > 20.96) return .58417 * pow(att -20.96, .4) + .07886 * (att - 20.96);
  return 0;
}

static
double * lsx_make_lpf(int num_taps, double Fc, double beta, double rho,
    double scale, sox_bool dc_norm)
{
  int i, m = num_taps - 1;
  double * h = malloc(num_taps * sizeof(*h)), sum = 0;
  double mult = scale / lsx_bessel_I_0(beta), mult1 = 1 / (.5 * m + rho);
  assert(Fc >= 0 && Fc <= 1);

  for (i = 0; i <= m / 2; ++i) {
    double z = i - .5 * m, x = z * M_PI, y = z * mult1;
    h[i] = x? sin(Fc * x) / x : Fc;
    sum += h[i] *= lsx_bessel_I_0(beta * sqrt(1 - y * y)) * mult;
    if (m - i != i)
      sum += h[m - i] = h[i];
  }
  for (i = 0; dc_norm && i < num_taps; ++i) h[i] *= scale / sum;
  return h;
}

static
void lsx_kaiser_params(double att, double Fc, double tr_bw, double * beta, int * num_taps)
{
  *beta = *beta < 0? lsx_kaiser_beta(att, tr_bw * .5 / Fc): *beta;
  att = att < 60? (att - 7.95) / (2.285 * M_PI * 2) :
    ((.0007528358-1.577737e-05**beta)**beta+.6248022)**beta+.06186902;
  *num_taps = !*num_taps? ceil(att/tr_bw + 1) : *num_taps;
}

double * lsx_design_lpf(
    double Fp,      /* End of pass-band */
    double Fs,      /* Start of stop-band */
    double Fn,      /* Nyquist freq; e.g. 0.5, 1, PI */
    double att,     /* Stop-band attenuation in dB */
    int * num_taps, /* 0: value will be estimated */
    int k,          /* >0: number of phases; <0: num_taps ß 1 (mod -k) */
    double beta)    /* <0: value will be estimated */
{
  int n = *num_taps, phases = max(k, 1), modulo = max(-k, 1);
  double tr_bw, Fc, rho = phases == 1? .5 : att < 120? .63 : .75;

  Fp /= fabs(Fn), Fs /= fabs(Fn);        /* Normalise to Fn = 1 */
  tr_bw = .5 * (Fs - Fp); /* Transition band-width: 6dB to stop points */
  tr_bw /= phases, Fs /= phases;
  tr_bw = min(tr_bw, .5 * Fs);
  Fc = Fs - tr_bw;
  assert(Fc - tr_bw >= 0);
  lsx_kaiser_params(att, Fc, tr_bw, &beta, num_taps);
  if (!n)
    *num_taps = phases > 1 ? *num_taps / phases * phases + phases - 1
                           : (*num_taps + modulo - 2) / modulo * modulo + 1;
  return Fn < 0? 0 : lsx_make_lpf(
      *num_taps, Fc, beta, rho, (double)phases, sox_false);
}
