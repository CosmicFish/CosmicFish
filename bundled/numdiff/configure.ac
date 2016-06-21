dnl
dnl Autoconf configuration script for Numdiff - version 5.x  .
dnl
dnl Copyright 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013  Ivano Primi  <ivprimi(at)libero(dot)it>.
dnl
dnl Permission is granted to use, copy, modify, and redistribute this file as
dnl needed.
dnl
dnl Process this file with autoconf to produce a configure script.
dnl

dnl --------------------------------------
dnl Minimal required version of 
dnl autoconf
dnl --------------------------------------
dnl

AC_PREREQ(2.62)

dnl --------------------------------------
dnl Initialization
dnl --------------------------------------
dnl

AC_INIT(Numdiff, 5.8, <ivprimi at libero dot it>, numdiff)
AC_CONFIG_SRCDIR(numdiff.h)
AC_CONFIG_HEADER(config.h)

dnl --------------------------------------
dnl My own flags
dnl --------------------------------------
dnl
if test "x$DBGFLAGS" = "x"; then
   DBGFLAGS="-g"
fi

if test "x$OPTFLAGS" = "x"; then
   OPTFLAGS="-O"
fi

dnl --------------------------------------
dnl Getting environment variables
dnl --------------------------------------
dnl
CPPFLAGS="${CPPFLAGS:=}"
CFLAGS="${CFLAGS:=}"
LDFLAGS="${LDFLAGS:=}"
LIBS="${LIBS:=}"
CPPFLAGS="$CPPFLAGS -DHAVE_CONFIG_H"

dnl --------------------------------------
dnl Check whether to add /usr/local or not
dnl (this is somewhat a religious problem)
dnl --------------------------------------
dnl
if test "`$CPP -v < /dev/null 2>&1 | grep '/usr/local/include' 2>&1`" = ""; then
  CPPFLAGS="$CPPFLAGS -I/usr/local/include"
  LDFLAGS="$LDFLAGS -L/usr/local/lib"
fi

dnl --------------------------------------
dnl If possible, enable
dnl extensions to C or Posix on hosts 
dnl that normally disable the extensions,
dnl typically due to standards-conformance
dnl namespace issues.  This should be 
dnl called before any macros that run the
dnl C compiler.
dnl --------------------------------------

AC_USE_SYSTEM_EXTENSIONS

dnl --------------------------------------
dnl Process options
dnl --------------------------------------
dnl
AC_ARG_ENABLE(debug, [  --enable-debug          turn on debugging  [default=no]],
                   , enable_debug=no)

AC_ARG_ENABLE(optimization, [  --enable-optimization   turn on basic optimization  [default=no]],
                          , enable_optimization=no)

AC_ARG_ENABLE(nls, [  --enable-nls            turn on Natural Language Support  [default=yes]],
                 , enable_nls=yes)

AC_ARG_ENABLE(gmp, [  --enable-gmp            turn on use of GNU MP  [default=yes]],
                 , enable_gmp=yes)

dnl --------------------------------------
dnl Checks for programs
dnl --------------------------------------
dnl
AC_PROG_CC
AC_PROG_CC_C_O
AC_PROG_CPP
dnl AC_PROG_INSTALL
AC_PATH_PROG(FIND,find)
dnl * Debian uses ginstall-info instead of install-info
dnl * to install info files.
AC_PATH_PROGS(INSTINFO,install-info ginstall-info)

dnl --------------------------------------
dnl Checks for header files
dnl --------------------------------------
dnl
AC_HEADER_STDC
if test $ac_cv_header_stdc = yes; then
 AC_MSG_RESULT(Ok! I have just found them.)
else
  AC_MSG_ERROR(Cannot find ANSI C header files!)
fi

AC_HEADER_STDBOOL
AC_CHECK_HEADERS(fcntl.h locale.h stddef.h \
  stdlib.h string.h sys/file.h time.h unistd.h wchar.h wctype.h)
AC_HEADER_DIRENT
AC_HEADER_STAT
AC_HEADER_SYS_WAIT

dnl --------------------------------------
dnl Checks for typedefs, structures, and 
dnl compiler characteristics
dnl --------------------------------------
dnl
dnl AC__GNU_SOURCE
dnl AC_DEFINE([__EXTENSIONS__], 1, [Enable Solaris extensions.])
AC_SYS_LARGEFILE

AC_C_CONST
AC_C_INLINE
AC_C_VARARRAYS

AC_TYPE_SIZE_T
AC_TYPE_PID_T
AC_CHECK_TYPE(ptrdiff_t, int)
AC_CHECK_TYPE(ssize_t, int)
AC_STRUCT_ST_BLKSIZE
AC_STRUCT_ST_RDEV
AC_STRUCT_ST_MTIM_NSEC
dnl AC_TYPE_MODE_T
dnl AC_STRUCT_ST_BLOCKS

dnl --------------------------------------
dnl Checks for libraries
dnl --------------------------------------
if test "x$enable_gmp" = "xyes"; then
  AC_CHECK_LIB(gmp, __gmpz_init)
fi

dnl Required to make building under Cygwin possible
if test "x$enable_nls" = "xyes"; then
  AC_CHECK_LIB(intl, setlocale)
fi

dnl --------------------------------------
dnl Checks for library functions
dnl --------------------------------------
dnl
dnl HAVE_SNPRINTF=0
dnl AC_CHECK_FUNC(snprintf, HAVE_SNPRINTF=1)
dnl if test $HAVE_SNPRINTF = 0; then
dnl   AC_MSG_ERROR(Declaration of snprintf() not present!)
dnl else
dnl   AC_MSG_RESULT(Ok! snprintf() is declared.)
dnl fi

HAVE_VSPRINTF=0
AC_CHECK_FUNC(vsprintf, HAVE_VSPRINTF=1)
if test $HAVE_VSPRINTF = 0; then
  AC_MSG_ERROR(Declaration of vsprintf() not present!)
else
  AC_MSG_RESULT(Ok! vsprintf() is declared.)
fi

HAVE_MEMSET=0
AC_CHECK_FUNC(memset, HAVE_MEMSET=1)
if test $HAVE_MEMSET = 0; then
  AC_MSG_ERROR(Declaration of memset() not present!)
else
  AC_MSG_RESULT(Ok! memset() is declared.)
fi

HAVE_MEMCHR=0
AC_CHECK_FUNC(memchr, HAVE_MEMCHR=1)
if test $HAVE_MEMCHR = 0; then
  AC_MSG_ERROR(Declaration of memchr() not present!)
else
  AC_MSG_RESULT(Ok! memchr() is declared.)
fi

HAVE_STRERROR=0
AC_CHECK_FUNC(strerror, HAVE_STRERROR=1)
if test $HAVE_STRERROR = 0; then
  AC_MSG_ERROR(Declaration of strerror() not present!)
else
  AC_MSG_RESULT(Ok! strerror() is declared.)
fi

HAVE_STRSPN=0
AC_CHECK_FUNC(strspn, HAVE_STRSPN=1)
if test $HAVE_STRSPN = 0; then
  AC_MSG_ERROR(Declaration of strspn() not present!)
else
  AC_MSG_RESULT(Ok! strspn() is declared.)
fi

HAVE_STRCSPN=0
AC_CHECK_FUNC(strcspn, HAVE_STRCSPN=1)
if test $HAVE_STRCSPN = 0; then
  AC_MSG_ERROR(Declaration of strcspn() not present!)
else
  AC_MSG_RESULT(Ok! strcspn() is declared.)
fi

AC_FUNC_STRERROR_R
AC_FUNC_FORK
AC_CHECK_FUNCS(localeconv gettimeofday isascii setlocale strchr strftime waitpid)

# setmode DOS ?
AC_FUNC_SETMODE_DOS

# Check for clock_gettime and its library.  Solaris puts it in -lrt or
# -lposix4, but we don't want to link that library unless we have to.
diff_saved_libs=$LIBS
AC_SEARCH_LIBS(clock_gettime, [rt posix4])
case $ac_cv_search_clock_gettime in
no | 'none required')
  LIB_CLOCK_GETTIME=;;
*)
  LIB_CLOCK_GETTIME=$ac_cv_search_clock_gettime;;
esac
AC_SUBST([LIB_CLOCK_GETTIME])
AC_CHECK_FUNCS(clock_gettime)
LIBS=$diff_saved_libs

dnl --------------------------------------
dnl Check for functions declarations
dnl --------------------------------------
dnl
AC_DEFINE([USE_UNLOCKED_IO], [1],
     [Define to 1 if you want getc etc. to use unlocked I/O if available.
      Unlocked I/O can improve performance in unithreaded apps,
      but it is not safe for multithreaded apps.])

AC_CHECK_DECLS(
     [clearerr_unlocked, feof_unlocked, ferror_unlocked,
      fflush_unlocked, fgets_unlocked, fputc_unlocked, fputs_unlocked,
      fread_unlocked, fwrite_unlocked, getc_unlocked,
      getchar_unlocked, putc_unlocked, putchar_unlocked, strerror])

dnl --------------------------------------
dnl Optimization, debug, GNU MP and NLS
dnl --------------------------------------
dnl

if test "x$enable_gmp" = "xno"; then
  CPPFLAGS="$CPPFLAGS -DDISABLE_GMP "
fi

if test "x$enable_nls" = "xyes"; then
  CPPFLAGS="$CPPFLAGS -DENABLE_NLS "
fi

if test "x$enable_debug" = "xyes"; then
   CFLAGS="$CFLAGS $DBGFLAGS "
fi

if test "x$enable_optimization" = "xyes"; then
   CFLAGS="$CFLAGS $OPTFLAGS "
fi

dnl LIBS="$LIBS"
dnl CPPFLAGS="$CPPFLAGS"

dnl --------------------------------------
dnl Output the makefile
dnl --------------------------------------
dnl
AC_MSG_RESULT(Configuration terminated with success:)
AC_MSG_RESULT(Preprocessor flags = "$CPPFLAGS")
AC_MSG_RESULT(Compiler flags     = "$CFLAGS")
AC_MSG_RESULT(Linker flags       = "$LDFLAGS")
AC_MSG_RESULT(Libraries          = "$LIBS")
dnl Remember to include @SET_MAKE@ in Makefile.in !!!
AC_PROG_MAKE_SET
AC_OUTPUT(Makefile)

dnl
dnl --------------------------------------
dnl End of configure script.
dnl --------------------------------------
