// -*- C++ -*-
// CLASSDOC OFF
// $Id: Transformation.h 478 2010-01-22 08:54:39Z stroili $
// ---------------------------------------------------------------------------
// CLASSDOC ON
//
// This file is part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the definition of the HepTransformation class for performing
// transformations on objects of the Hep3Vector and HepPoint
// classes
//
// .SS See Also
// ThreeVector.h, HepPoint.h, Rotation.h, Translation.h
//
// .SS History
// Author Victoria Novotny (LBL)
// Modified by Gautier Hamel de Monchenault & David N Brown

#ifndef _TRANSFORMATION_H_
#define _TRANSFORMATION_H_

class HepPoint;
#include "CLHEP/Vector/ThreeVector.h"
#include "BTrk/BbrGeom/Translation.h"
#include "CLHEP/Vector/Rotation.h"
class EulerAngles;
class AlignAngles;
using namespace CLHEP;

//
// A transformation takes a rigid body, translates it
// by a translation vector trans_t and rotates it by a rotation matrix trans_r.
//
class HepTransformation{

public:
  ~HepTransformation();
// destructor

  HepTransformation();
// default constructor

  HepTransformation( const HepRotation&, const HepTranslation& );
// constructor

  HepTransformation( const Hep3Vector&, const Hep3Vector& );
// vector and normal

  HepTransformation( const Hep3Vector&, const Hep3Vector&, const Hep3Vector&);
// vector, normal and x axis 

  HepTransformation( const Hep3Vector&,const EulerAngles&);
// vector and 3 Euler angles

  HepTransformation( const Hep3Vector&,const AlignAngles&);
// vector and 3 Alignment (small) angles

  HepTransformation(const HepTransformation &);
// copy constructor

  HepTransformation & operator = (const HepTransformation & t);
// Assignment

  HepTransformation& operator *= (const HepTransformation&);
  HepTransformation& transform(const HepTransformation&);
// concatonation operators.  *= describes the net transformation
// applying the second to the origin of the first.  'transform'
// describes the reverse order.
  HepTransformation& invert(); // invert the transform in place
  HepTransformation inverse(); // return the inverse transform

  HepRotation rot_mat() const;  // returns Rotation matrix

  HepTranslation tra_vec()  const;  // returns Translation vector
//
// coordinate transforms
//
  HepPoint    transFrom( const HepPoint&   ) const;
  HepPoint    transTo  ( const HepPoint&   ) const;
  Hep3Vector  transFrom( const Hep3Vector& ) const;
  Hep3Vector  transTo  ( const Hep3Vector& ) const;
// origin and basis vector of the second frame
//  (these ones won't be used very often - no need to return a reference)
  HepPoint   origin() const;       
  Hep3Vector unit( int axis=2 ) const; 

// Transform vector "in place"
  Hep3Vector & transform(Hep3Vector &) const;

private:

  HepTranslation trans_t;
// Translation vector components

  HepRotation trans_r;
// Rotation matrix components

};

#endif

