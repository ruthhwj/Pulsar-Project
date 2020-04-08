noinst_LTLIBRARIES=libsalsa.la

libsalsa_la_SOURCES=psrsalsa_defines.h psrsalsa.h psrsalsa_typedefs.h amoeba_d.c amoeba_nmsimplex_d.c angles.c application.c astronomy.c fft.c fitting.c fold.c linalg.c minimize.c myio.c pgplot.c ppgplot.c psrio.c psrio_epn.c psrio_fits.c psrio_opperations.c psrio_paswing.c psrio_preprocess.c psrio_puma.c psrio_sigproc.c pulseenergy.c rmsynth.c spectra.c statistics.c vonMises.c psrsalsa_tempo3.h amoeba.c amoeba_ld.c amoeba_nmsimplex.c amoeba_nmsimplex_d.c amoeba_nmsimplex_ld.c binarymodels.c binarymodels_simplyfied.c calc.c fitting.c fitting_ld.c linalg_ld.c pgplotsubroutineshelp.c psrio_desh.c psrio_fb.c psrio_gmrt.c psrio_presto.c tempo3_lib.c cal.c mpolyco.c ppolyco.c

libsalsa_la_LIBADD=../slalib/libslawrap.la @FFTW_LIBS@ @CFITSIO_LIBS@ @PGPLOT_LIBS@ @GSL_LIBS@

AM_CPPFLAGS=-I$(top_srcdir)/src/slalib @FFTW_CFLAGS@ @CFITSIO_CFLAGS@ @PGPLOT_CFLAGS@ @GSL_CFLAGS@ @GTK_CFLAGS@

#if NEW_GSL

#AM_CPPFLAGS+=-DGSL_VERSION_2_1

#endif

