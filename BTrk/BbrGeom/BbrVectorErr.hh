//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: BbrVectorErr.hh 491 2010-01-13 16:59:16Z stroili $
//
// Description:
//	Add errors to a vector.  Used for direction errors
//      BaBar native class
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Forest Rouse	                February 1996
//      Victoria Novotny                August   1996
//      Ed Frank        University of Pennsylvania, efrank@upenn5.hep.upenn.edu
//
// History
//      14 Oct 96  Ed Frank      Simple mod to make unary - compile under
//                               Sun's CC 4.0.1
//      19 Jan 2002  Sasha Telnov  Added operator * (scaling by a real number)
//
// Copyright Information:
//	Copyright (C) 1996
//
//------------------------------------------------------------------------
#ifndef BBRVECTORERR_HH
#define BBRVECTORERR_HH

#include <iosfwd>
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BbrGeom/BbrError.hh"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
// #include "BTrk/BbrGeom/Translation.h"
// #include "BTrk/BbrGeom/Transformation.h"
using namespace CLHEP;

class BbrVectorErr : public Hep3Vector {
   public:
    // polar coordinates
    enum PolarCoordinateIndex { Rho = 0, Theta = 1, Phi = 2, NUM_PCOORDINATES = 3 };

    enum CylindricalCoordinateIndex { C_Rho = 0, C_Zeta = 1, C_Phi = 2, NUM_CCOORDINATES = 3 };
    // argumentless constructor:
    BbrVectorErr() : Hep3Vector(), _covMatrix(NUM_COORDINATES) {}

    // auto casting constructor
    BbrVectorErr(const Hep3Vector& p) : Hep3Vector(p), _covMatrix(NUM_COORDINATES) {}
    BbrVectorErr(const Hep3Vector& p, const BbrError& covMat)
        : Hep3Vector(p), _covMatrix(NUM_COORDINATES) {
        _covMatrix = covMat;
    }

    // copy constructor:
    BbrVectorErr(const BbrVectorErr& v) : Hep3Vector(v), _covMatrix(v.covMatrix()) {}

    // destructor MAY be needed later
    // virtual ~BbrVectorErr() {};

    // assignment operator:
    BbrVectorErr& operator=(const BbrVectorErr& v) {
        if (this != &v) {
            Hep3Vector::operator=(v);
            _covMatrix = v.covMatrix();
        }
        return *this;
    }

    BbrVectorErr operator-() {
        Hep3Vector t = *this;
        return BbrVectorErr(-t, _covMatrix);  // _covMatrix remains unaltered
    }

    BbrVectorErr& operator+=(const BbrVectorErr& v) {
        Hep3Vector::operator+=(v);
        _covMatrix += v.covMatrix();
        return *this;
    }

    BbrVectorErr& operator-=(const BbrVectorErr& v) {
        Hep3Vector::operator-=(v);
        _covMatrix += v.covMatrix();
        return *this;
    }

    BbrVectorErr& transform(const HepTranslation& trans) {
        trans.transform(*this);
        return *this;
    }

    BbrVectorErr& transform(const HepRotation& rot) {
        Hep3Vector::transform(rot);
        _covMatrix = _covMatrix.similarity(rot);
        return *this;
    }

    BbrVectorErr& transform(const HepTransformation& transf) {
        transf.transform(*this);
        _covMatrix = _covMatrix.similarity(transf.rot_mat());
        return *this;
    }

    double determineChisq(const Hep3Vector& refVector) const;
    // returns Chisquare
    // refVector refers to the same origin as the Hep3Vector of this
    // ie refVector is not relative to this Vector

    inline const BbrError& covMatrix() const { return _covMatrix; }

    BbrError covRTPMatrix() const;
    // returns the covariance Matrix in spherical coordinate
    // use   PolarCoordinateIndex enum to get the components

    BbrError covRZPMatrix() const;
    // returns the covariance Matrix in cylindrical coordinate
    // use   CylindricalCoordinateIndex enum to get the components

    inline void setCovMatrix(const BbrError& v) { _covMatrix = v; }

    //  void printOn(ostream& out=cout) const;

   private:
    BbrError _covMatrix;
};

BbrVectorErr operator+(const BbrVectorErr&, const BbrVectorErr&);

BbrVectorErr operator-(const BbrVectorErr&, const BbrVectorErr&);

// Added by Sasha Telnov
BbrVectorErr operator*(const BbrVectorErr&, HepDouble a);
BbrVectorErr operator*(HepDouble a, const BbrVectorErr&);

std::ostream& operator<<(std::ostream& stream, const BbrVectorErr& verr);
std::istream& operator>>(std::istream& stream, BbrVectorErr& verr);

#endif
