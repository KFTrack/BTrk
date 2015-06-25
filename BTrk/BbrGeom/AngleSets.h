// -*- C++ -*-
// CLASSDOC OFF
// $Id: AngleSets.h 478 2010-01-22 08:54:39Z stroili $
// ---------------------------------------------------------------------------
// CLASSDOC ON
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

#ifndef ANGLESETS_H
#define ANGLESETS_H

class EulerAngles{
public:
  double a,b,c; // 3 euler angles
  EulerAngles(double,double,double);
};

class AlignAngles{
public:
  double a,b,c; // 3 angles for succesive rotations about the x,y, and z axes
  AlignAngles(double,double,double);
};

#endif
