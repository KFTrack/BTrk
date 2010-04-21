// -*- C++ -*-
// CLASSDOC OFF
// $Id: Translation.h 478 2010-01-22 08:54:39Z stroili $
// ---------------------------------------------------------------------------
// CLASSDOC ON
//
// This file is part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the definition of the HepTranslation class for performing 
// translations on objects of the Hep3Vector and HepPoint classes
//
// .SS See Also
// ThreeVector.h, HepPoint.h
//
// .SS History
// Author Victoria Novotny (LBL)

#ifndef _TRANSLATION_H_
#define _TRANSLATION_H_

#ifdef __GNUC__
#pragma interface
#endif

#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/HepPoint.h"

#ifdef HEP_NO_INLINE_IN_DECLARATION
#define inline
#endif

class HepPoint;
#include "CLHEP/Vector/ThreeVector.h"

class HepTranslation{

public:
    ~HepTranslation();
    // the destructor
    
    HepTranslation();
    // default constructor

    HepTranslation(const Hep3Vector &);
    // constructor

    HepTranslation(const HepTranslation &);
    // copy constructor

    HepTranslation & operator = (const HepTranslation &);
    // Assignment

    double x() const;
    double y() const;
    double z() const;
    // Get the x, y, z components
    
    Hep3Vector trans_vec() const;
    // Returns the vector components

    HepTranslation inverse() const;
    // Returns inverse

    HepTranslation & invert();
    // Inverts the Translation vector

    HepTranslation & translateX(double);
    // Translate along x-axis

    HepTranslation & translateY(double);
    // Translate along y-axis

    HepTranslation & translateZ(double);
    // Translate along z-axis
    
    HepTranslation & operator += (const HepTranslation &);
    // Addition of two translations

    Hep3Vector & transform(Hep3Vector &) const;
    // Transform vector "in place"

private:

    Hep3Vector vec;
    // Vector containing components

};

HepTranslation operator + (const HepTranslation &, const HepTranslation &);
// Addition Translation + Translation




#endif

