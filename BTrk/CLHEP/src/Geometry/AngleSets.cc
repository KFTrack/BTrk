// -*- C++ -*-
// $Id: AngleSets.cc 192 2009-03-04 12:20:53Z stroili $
// ---------------------------------------------------------------------------
//
// This file is part of the CLHEP - a Class Library for High Energy Physics.
//
// Classes for describing the two types of angles.  As their representation is
// identical (3 doubles), they must be distinguished as classes and not just
// typedefs in order to differentiate the respective constructors which use
// them.
//
// History:
// 02.07.97 A.Mokhtarani - initial version

#include "CLHEP/Geometry/AngleSets.h"

EulerAngles::EulerAngles(double ea,double eb,double ec)
  : a(ea),b(eb),c(ec)
{;}

AlignAngles::AlignAngles(double ea,double eb,double ec)
  : a(ea),b(eb),c(ec)
{;}
