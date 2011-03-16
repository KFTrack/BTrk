// -*- C++ -*-
// CLASSDOC OFF
// $Id: TemplateFunctions.h 192 2009-03-04 12:20:53Z stroili $
// ---------------------------------------------------------------------------
// CLASSDOC ON
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This file contains template definitions of some usefull functions.
//
#ifndef HEP_TEMPLATEFUNCTIONS_H
#define HEP_TEMPLATEFUNCTIONS_H

#include "CLHEP/config/CLHEP.h"
#include <algorithm>
using std::min;
using std::max;

#ifndef CLHEP_SQR_DEFINED
#define CLHEP_SQR_DEFINED
#ifdef sqr
#undef sqr
#endif

#ifndef HEP_SQR_NEEDS_PARAMETER_WITHOUT_CONST
template <class T>
inline T sqr(const T& x) {
  return x*x;
}
#else
template <class T>
inline T sqr(T x) {
  return x*x;
}
#endif /* HEP_SQR_NEEDS_PARAMETER_WITHOUT_CONST */
#endif

#ifndef CLHEP_ABS_DEFINED
#define CLHEP_ABS_DEFINED
#ifdef abs
#undef abs
#endif

#ifndef HEP_ABS_NEEDS_PARAMETER_WITHOUT_CONST
template <class T>
inline T abs(const T& a) {
  return a < 0 ? -a : a;
}
#else
template <class T>
inline T abs(T a) {
  return a < 0 ? -a : a;
}
#endif /* HEP_ABS_NEEDS_PARAMETER_WITHOUT_CONST */
#endif

#endif /* HEP_TEMPLATEFUNCTIONS_H */








