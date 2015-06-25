// -*- C++ -*-
// CLASSDOC OFF
// $Id: HepPoint.h 478 2010-01-22 08:54:39Z stroili $
// ---------------------------------------------------------------------------
// CLASSDOC ON
//
// This file is part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the definition of the HepPoint class.
//
// .SS See Also
// ThreeVector.h, Rotation.h, Translation.h
//
// .SS History
// Author:Victoria Novotny(LBL)

#ifndef _HEPPOINT_H_
#define _HEPPOINT_H_

#include "CLHEP/Vector/Rotation.h"
class HepTranslation;
class HepTransformation;
#include "CLHEP/Vector/ThreeVector.h"
class HepPoint{

public:
  enum CoordinateIndex {   
    X = 0, 
    Y = 1, 
    Z = 2,
    NUM_COORDINATES = 3,
    SIZE = NUM_COORDINATES
  };
  // Safe indexing of the coordinates when using with matrices, arrays, etc.

    HepPoint(double x = 0.0, double y = 0.0, double z = 0.0);
    // Constructor

    HepPoint(const HepPoint &);
    // Constructor

    inline ~HepPoint(){}
    // Destructor

    HepPoint & operator = (const HepPoint &);
    // Assignment

    HepPoint operator - () const;
    //Unary minus

    HepPoint & operator += (const CLHEP::Hep3Vector &);
    // Addition

    HepPoint & operator -= (const CLHEP::Hep3Vector &);
    // Subtraction

    inline double x() const{ return dx; }
    inline double y() const{ return dy; }
    inline double z() const{ return dz; }
    // Get the x, y, z components

    double mag2() const;
    // The magnitude squared

    double mag() const;
    // The magnitude

    double phi() const;
    // The azimuth angle

    double theta() const;
    // The polar angle

    double cosTheta() const;
    // Cosine of polar angle
  
    double distanceTo(const HepPoint&) const;
    // distance from this point to another point
  
    HepPoint & operator *= (const CLHEP::HepRotation &);
    HepPoint & transform(const CLHEP::HepRotation &);
    // Transformation with a Rotation matrix. Note a*=b; <=> a=b*a;

    HepPoint & transform(const HepTranslation &);
    // Transform a point with a Translation vector

    HepPoint & transform(const HepTransformation &);
    // Transform a point with a Transformation matrix
 
    void rotate(float angle, const CLHEP::Hep3Vector & axis);
    // Rotates HepPoint around axis specified by a CLHEP::Hep3Vector

    void rotateX(float);
    // Rotates HepPoint around the x-axis

    void rotateY(float);
    // Rotates HepPoint around the y-axis

    void rotateZ(float);
    // Rotates HepPoint around the z-axis

    inline void setX(double x){ dx = x; }
    inline void setY(double y){ dy = y; }
    inline void setZ(double z){ dz = z; }
    // Set the x, y, z components

    void setMag(double);
    void setTheta(double);
    void setPhi(double);
    // Change x, y, z components according to new magnitude, theta and phi 

private:

    double dx, dy, dz;
    // The x, y, z components

};


std::istream & operator >> (std::istream &, HepPoint &);
std::ostream & operator << (std::ostream &, const HepPoint &);
// output(intput) to(from) a stream


int operator == (const HepPoint &, const HepPoint &);
int operator != (const HepPoint &, const HepPoint &);
// Equality operators

HepPoint operator + (const HepPoint &, const CLHEP::Hep3Vector &);
// Addition point + vector

HepPoint operator + (const CLHEP::Hep3Vector &, const HepPoint &);
// Addition vector + point

CLHEP::Hep3Vector operator - (const HepPoint &, const HepPoint &);

HepPoint operator - (const HepPoint &, const CLHEP::Hep3Vector &);
// Subtraction point - vector

HepPoint operator * (const CLHEP::HepRotation &, const HepPoint &);

#endif

