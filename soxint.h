/* 
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
#ifndef SOXINT_H
#define SOXINT_H

#include <sys/types.h>
#include <stddef.h>
#include <stdlib.h> /* might have min() and max(), but not always */
#include <stdint.h>
#include <malloc.h>

#ifndef max
# define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
# define min(a,b) (((a) < (b)) ? (a) : (b))
#endif
#define range_limit(x, lower, upper) (min(max(x, lower), upper))

#define array_length(a) (sizeof(a)/sizeof(a[0]))

typedef enum sox_bool {
    sox_false,
    sox_true
} sox_bool;

int _lsx_set_dft_length(int num_taps);

#if defined(_MSC_VER) && _MSC_VER >= 1400
# define HAVE__ALIGNED_MALLOC 1
#endif

#if HAVE_MEMALIGN
# define lsx_aligned_malloc(align, size) memalign(align, size)
# define lsx_aligned_free free
#elif HAVE__ALIGNED_MALLOC
# define lsx_aligned_malloc(align, size) _aligned_malloc(size, align)
# define lsx_aligned_free _aligned_free
#else
# define lsx_aligned_malloc(align, size) pffft_aligned_malloc(size)
# define lsx_aligned_free pffft_aligned_free
#endif

#endif
