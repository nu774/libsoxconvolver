ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = subdir-objects

lib_LTLIBRARIES = libsoxconvolver.la
noinst_LTLIBRARIES = libpffft_ex.la

libsoxconvolver_la_SOURCES = \
    dft_filter.c    \
    effects_i_dsp.c \
    fft4g_single.c

libpffft_ex_la_SOURCES = pffft_ex.c

if HAVE_WINDRES
    libsoxconvolver_la_SOURCES += libsoxconvolver.rc
endif

libpffft_ex_la_CFLAGS = $(AM_CFLAGS) -prefer-pic

if GCC
if IX86
libpffft_ex_la_CFLAGS += -msse -mfpmath=sse
if HAVE_WINDRES
libpffft_ex_la_CFLAGS += -mincoming-stack-boundary=2
endif
endif
endif

libsoxconvolver_la_LDFLAGS = -no-undefined -export-symbols-regex '^lsx_.*$$' \
			     -version-info @VI_CURRENT@:@VI_REVISION@:@VI_AGE@
libsoxconvolver_la_LIBADD = libpffft_ex.la $(LIBM)

AM_CPPFLAGS = -DNDEBUG

include_HEADERS = libsoxconvolver.h

.rc.lo: 
	$(LIBTOOL) $(AM_LIBATOOLFLAGS) $(LIBTOOLFLAGS) \
	    --tag=RC --mode=compile $(RC) $(RCFLAGS) $< -o $@

