#ifndef BABAR_HH
#define BABAR_HH
//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: BaBar.hh 588 2010-02-25 12:04:05Z stroili $
//
// Description:
//	Standard include file to contain global "fixes"
//      and type definitions
//
//      This file is intended to be included by _lots_ of others
//      so should stay very light-weight.  In particular, definitions
//      of static variables are not appropriate here, as they
//      will allocate memory in every BaBar object file.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	David Quarrie
//      Bob Jacobsen (added F90 interface types, fixed bool, Mar 19,95)
//	Forest Rouse (added pi, bznom, chisq_undef, and epsilon)
//      Abi Soffer (added BbrLorentzVector)
//      A. De Silva (moved macros from Framework/AppModule.hh. Jan 29,99)
//
// Copyright Information:
//	Copyright (C) 1995, 1996
//
//------------------------------------------------------------------------

// Need to make sure the very minimum gets into the dataflow builds
#ifndef VXWORKS

// Define very large integer types. These types are required by the
// POSIX specification.

#if defined(__SUNPRO_CC) || defined(__linux__)
// ROOT's CINT cannot handle the type long long
#ifndef __CINT__
#if defined(__x86_64__)
#include <stdlib.h>
#include <stdint.h>
#else
typedef long long                 int64_t;
typedef unsigned long long       uint64_t;
#endif
#endif
#endif

// Activate large file support for Sun if it isn't already.
// This should be done in CPPFLAGS (and is) but a bug in the template
// support in the Workshop 4.2 compiler causes some files to be
// compiled with _FILE_OFFSET_BITS set to "'64'" instead of "64".  // FIXME
// NB: This _MUST_ be done before <sys/feature_tests.h> is included.
#if defined(__SUNPRO_CC)
#ifdef _LARGEFILE_SOURCE
#undef _LARGEFILE_SOURCE
#endif
#define _LARGEFILE_SOURCE 1
#ifdef _FILE_OFFSET_BITS
#undef _FILE_OFFSET_BITS
#endif
#define _FILE_OFFSET_BITS 64
#endif // Sun


// The following two "paragraphs" of code seem to make the ROOT
// interpreter unhappy.  For now we protect it from them.  This is
// not an ideal solution.                                        // FIXME
#ifndef __CINT__

// Load the various platforms' feature-test include files, so that we
// can properly take advantage of the various *POSIX* and *XOPEN* macros.
// There is no way to do this in a platform-independent way without 
// loading something heavyweight like <unistd.h>, which might affect
// code that doesn't want it.  [gpdf 2001.01.29]
#ifdef __SUNPRO_CC
 #include <sys/feature_tests.h>
#else
 #ifdef __linux__
 #include <features.h>
  #else
   // Don't know what to do here.  Try loading a lightweight include file.
   // On some platforms this brings in the feature tests, but not on all.
   #include <stddef.h>
  #endif
#endif


// Accommodate the variations in the adoption of socklen_t from Unix98.
// Unix98 (XPG5) has socklen_t.  Unix95 (XPG4v2) uses size_t (which is
// typically unsigned) in places where socklen_t appears in Unix98.
// Pre-standard APIs tend to have int instead.    [added gpdf 2001.01.29]
#if ( defined(_XOPEN_SOURCE) && ( _XOPEN_SOURCE - 0 ) >= 500 )
// Unix98 (XPG5): Solaris 7 with -D_XOPEN_SOURCE=500, Linux
 #define BABAR_socklen_t socklen_t
#else
 #ifdef _XOPEN_SOURCE_EXTENDED
  // Unix95 (XPG4v2): OSF, Solaris 2.6 with -D_XOPEN_SOURCE_EXTENDED
  #define BABAR_socklen_t size_t
 #else
  // pre-standard: Solaris 7 & 2.6 "plain vanilla"
  #define BABAR_socklen_t int
 #endif
#endif

#endif // ifndef __CINT__

#if defined(__SUNPRO_CC) && \
      ( !defined(__EXTENSIONS__) && \
	( defined(_POSIX_C_SOURCE) || defined(_XOPEN_SOURCE) ) )
// This one is in the libraries on Solaris but only has a prototype
// defined in <signal.h> if pure-Unix95/98 compilation is not in force.
// We just restore the prototype.
extern "C" char* strsignal( int sig );
#endif

// An error in the construction of the system and C++ header files on
// Solaris 8 / Workshop 6 Updates 1&2 leads to a conflict between the use
// of ::clock_t and std::clock_t when <string> is compiled under
// -D_XOPEN_SOURCE=500.  The following code ensures that ::clock_t is 
// always defined and thus allows <string> to compile.  
// This is just a workaround and should be monitored as compiler and
// operating system versions evolve.  [gpdf 2002.02.05]
#if defined(__SUNPRO_CC) && defined(_XOPEN_SOURCE) && ( _XOPEN_SOURCE - 0 == 500 )
#ifndef _CLOCK_T
#define _CLOCK_T
typedef long            clock_t; /* relative time in a specified resolution */
#endif  /* ifndef _CLOCK_T */
#endif // SUN and XOPENSOURCE=500


#if defined(__SUNPRO_CC)
// There are two forms of each of the std::count and std::distance
// functions: a pre-standard version that requires an extra argument
// to hold the result, and a newer version, which made it into the
// standard, that returns its result.  On Solaris, the STL is built
// with the _RWSTD_NO_CLASS_PARTIAL_SPEC flag set, which disables the
// new forms.  The following is an emulation of the standard versions.
// Although this is not complete, it enables us to compile (our
// current) standard code on Solaris. [bartoldu 2004.07.20]
//
// I included std::count_if here as well which suffers from the same problem.
// [narsky 2005.02.18]

#include <stddef.h>

namespace std
{
  template <class InputIterator, class T>
  ptrdiff_t
  count (InputIterator first, InputIterator last, const T& value);

  template <class InputIterator, class Predicate>
  ptrdiff_t
  count_if (InputIterator first, InputIterator last, Predicate pred);

  template <class InputIterator>
  ptrdiff_t
  distance (InputIterator first, InputIterator last);
}

#ifdef    BABAR_COMP_INST
#include  "BTrk/BaBar/BaBar.icc"
#endif  // BABAR_COMP_INST

#endif // SUNCC


// This block is for Mac OSX 10.3.5 with gcc 3.3
#ifdef __APPLE_CC__

// No realtime signal support
#ifndef SIGRTMAX
#define SIGRTMAX SIGUSR1
#endif

// "union sigval" not properly defined
#define sigval_int sival_int
#define sigval_ptr sival_ptr

// No clock_gettime in <time.h>
typedef int clockid_t;
#ifndef CLOCK_REALTIME
#define CLOCK_REALTIME 0

// No ENODATA, this is used in OepRemoteFramework
#define ENODATA 61

// No O_LARGEFILE as this is a 64bit platform
// Used in OepFramework and OlmMerger
#define O_LARGEFILE 0

#endif

// <cmath> undefines isnan() and friends without providing replacements
// Only isnan() and isfinite() are used in BaBar; leave the others alone  

#if !defined(isfinite)
extern "C" int isfinite(double);
#endif
#if !defined(isnan)
extern "C" int isnan(double);
#endif

// ROOT's CINT has a problem with timespec not being defined
#ifndef __CINT__
int clock_gettime ( clockid_t clock_id, struct timespec *tp);
#endif

#endif // __APPLE_CC__

#endif // VXWORKS

// Only below here should be seen by VXWORKS

typedef double HepDouble;
typedef int    HepInt;
typedef float  HepFloat;
typedef bool   HepBoolean;


// Used by BTrk/DetectorModel/src/DetMaterial.cc
inline double sqr( double a){
  return a*a;
}

// Electron mass in GeV/c^2.  Needed by BTrk/DetectorModel/src/DetMaterial.cc
// Replace with proper invocation from PDT.
//static const double electron_mass_c2 = 0.000510998910;

#include <stdio.h>
#include <algorithm>

#endif // BABAR_HH
