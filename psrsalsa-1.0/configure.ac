#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([psrsalsa],[1.0],[patrick.weltevrede@manchester.ac.uk],[],[https://github.com/weltevrede/psrsalsa/])

AC_CONFIG_MACRO_DIR([autoconf])
AC_CONFIG_AUX_DIR([config])

# Enable large-file support
AC_SYS_LARGEFILE

# Determine the machine endian
AC_C_BIGENDIAN
AC_LANG([C])

# Checks for compilers
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_CC
AC_PROG_CPP
AM_PROG_CC_C_O
AC_PROG_F77
AC_PROG_F77_C_O
AC_F77_LIBRARY_LDFLAGS
AC_PROG_FC


AC_PROG_INSTALL
AC_PROG_LIBTOOL

# Check whether the 'cos' function can be found in library 'm'
# (standard C math library).
# Math library is required for the tempo2 predictor library
AC_CHECK_LIB(m, cos)
save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -Wno-error"
ACX_PTHREAD
CFLAGS="$save_CFLAGS"

AC_SUBST(PTHREAD_LIBS)
AC_SUBST(PTHREAD_CFLAGS)

MJK_FORTRAN_TEST

if test X$fortran_c_links != Xyes ; then
    AC_MSG_ERROR(Your C and F77 compilers do not cross link (\$CC=$CC \$F77=$F77)) 
fi

SWIN_LIB_PGPLOT
MJK_LIB_FFTW
SWIN_LIB_CFITSIO

#GSL basic stuff

AX_PATH_GSL(1.0.0, have_gsl=yes,have_gsl=no)
AM_CONDITIONAL(HAVE_GSL, [test x"$have_gsl" = xyes])
#gsl_version_integer = 100 * $gsl_major_version + $gsl_minor_version
#gsl_version_integer=`expr 100 \* $gsl_major_version + $gsl_minor_version`
gsl_version_integer=$((100 * $gsl_major_version + $gsl_minor_version))
AC_MSG_NOTICE(Found GSL version: $gsl_major_version.$gsl_minor_version = $gsl_version_integer)
GSL_CFLAGS="$GSL_CFLAGS -DGSL_VERSION_NUMBER=$gsl_version_integer"
AC_MSG_NOTICE(Found GSL CFLAGS: $GSL_CFLAGS)
#gsl_minor_version

#AC_CHECK_LIB(gsl,gsl_multifit_fdfsolver_jac,new_gsl=yes,new_gsl=no,"$GSL_LIBS")

#AM_CONDITIONAL(NEW_GSL, [test x"$new_gsl" = xyes])


#AM_PATH_GTK_2_0(2.24.0,have_gtk=yes,have_gtk=no)
AM_PATH_GTK_3_0(3.10.0,have_gtk=yes,have_gtk=no)
if test X$have_gtk != Xyes; then
    GTK_CFLAGS=""
else
    GTK_CFLAGS="$GTK_CFLAGS -DTEMPO3_EnableGTK=1"
fi
AC_MSG_NOTICE(Found GTK CFLAGS: $GTK_CFLAGS)

AC_F77_WRAPPERS

#ACX_BLAS(have_blas=yes,have_blas=no)
#AC_MSG_CHECKING([for BLAS])
#AC_MSG_RESULT($have_blas)
#AM_CONDITIONAL(HAVE_BLAS, [test x$have_blas = xyes])

# test for LAPACK acceleration
#ACX_LAPACK(have_lapack=yes, have_lapack=no)
#AC_MSG_CHECKING([for LAPACK])
#AC_MSG_RESULT($have_lapack)
#AM_CONDITIONAL(HAVE_LAPACK, [test x"$have_lapack" = xyes])

AM_INIT_AUTOMAKE([foreign subdir-objects])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/lib/Makefile
                 src/prog/Makefile
                 src/slalib/Makefile])

AC_OUTPUT

AC_MSG_NOTICE(Fortran compiler is $F77 $FFLAGS)
AC_MSG_NOTICE(C compiler is $CC $CFLAGS)
AC_MSG_NOTICE(Global linking flags are $LDFLAGS)


# User warnings
if test X$have_pgplot != Xyes; then
    AC_MSG_ERROR(PGPLOT was not detected) 
else
    AC_MSG_NOTICE(PGPLOT compile flags are "$PGPLOT_CFLAGS") 
    AC_MSG_NOTICE(PGPLOT linking flags are "$PGPLOT_LIBS") 
fi

if test X$have_gsl != Xyes; then
    AC_MSG_ERROR(GSL was not found) 
fi

if test X$have_fftw != Xyes; then
    AC_MSG_ERROR(FFTW was not found) 
fi

if test X$have_cfitsio != Xyes; then
    AC_MSG_ERROR(cfitsio was not found) 
fi
if test X$have_gtk != Xyes; then
    AC_MSG_NOTICE(gtk was not found) 
    AC_MSG_NOTICE(***************************************************************************)
    AC_MSG_NOTICE(WARNING: GTK3+ WAS NOT FOUND. PLEASE INSTALL GTK3+. OTHERWISE YOU COULD)
    AC_MSG_NOTICE(PROCEED WITH RUNNING MAKE ALTHOUGH THE RELATED FEATURES WILL BE DISABLED.)
    AC_MSG_NOTICE(***************************************************************************)
fi



#echo ""
#echo "****"
#echo "Ok! Good to start compiling... try make && make install'"
#echo "****"



