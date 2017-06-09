// -*- C++ -*-
// $Id: Transformation.cc 192 2009-03-04 12:20:53Z stroili $
// ---------------------------------------------------------------------------
//
// This file is part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the HepTransformation class.
//
// .SS History
// Author Victoria Novotny (LBL)
// Modified by Gautier Hamel de Monchenault & David N Brown

#ifdef __GNUC__
#pragma interface
#endif

#include <math.h>

#include "BTrk/BbrGeom/AngleSets.h"
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/BbrGeom/Transformation.h"
#include "CLHEP/Vector/RotationInterfaces.h"
#include "CLHEP/Vector/ThreeVector.h"

HepTransformation::~HepTransformation() {}

HepTransformation::HepTransformation() {}

HepTransformation::HepTransformation(const HepRotation& m, const HepTranslation& p) {
    trans_r = m;
    trans_t = p;
}

HepTransformation::HepTransformation(const HepTransformation& q) {
    trans_r = q.rot_mat();
    trans_t = q.tra_vec();
}

HepTransformation& HepTransformation::operator=(const HepTransformation& q) {
    if (this == &q)
        return *this;
    trans_r = q.rot_mat();
    trans_t = q.tra_vec();
    return *this;
}

HepTransformation::HepTransformation(const Hep3Vector& o, const Hep3Vector& z) : trans_t(o) {
    double ct = cos(z.theta());
    double st = sin(z.theta());
    double cp = cos(z.phi());
    double sp = sin(z.phi());
    trans_r = HepRotation(HepRep3x3(-ct * cp, sp, st * cp, -ct * sp, -cp, st * sp, st, 0., ct));
}

HepTransformation::HepTransformation(const Hep3Vector& o, const Hep3Vector& z, const Hep3Vector& x)
    : trans_t(o) {
    Hep3Vector k = z.unit();
    Hep3Vector i = (x - (x * k) * k).unit();
    Hep3Vector j = k.cross(i);
    trans_r = HepRotation(HepRep3x3(i.x(), j.x(), k.x(), i.y(), j.y(), k.y(), i.z(), j.z(), k.z()));
}

HepTransformation::HepTransformation(const Hep3Vector& o, const EulerAngles& eang)
    : trans_t(o), trans_r() {
    trans_r.rotateZ(eang.c);
    trans_r.rotateX(eang.b);
    trans_r.rotateZ(eang.a);
}

HepTransformation::HepTransformation(const Hep3Vector& o, const AlignAngles& aang)
    : trans_t(o), trans_r() {
    trans_r.rotateX(aang.a);
    trans_r.rotateY(aang.b);
    trans_r.rotateZ(aang.c);
}

HepTransformation& HepTransformation::operator*=(const HepTransformation& trans2) {
    //
    //  Rotate the 2nd transform translation vector back into the
    //  primary frame before adding translations.
    //
    HepTranslation rottrans(trans_r * trans2.trans_t.trans_vec());
    trans_t += rottrans;
    trans_r *= trans2.trans_r;
    return *this;
}

HepTransformation& HepTransformation::transform(const HepTransformation& trans2) {
    //
    //  Equivalent to the above, just reversing the roles of the two
    //  transforms.
    //
    HepTranslation rottrans(trans2.trans_r * trans_t.trans_vec());
    trans_t = rottrans + trans2.trans_t;
    trans_r.transform(trans2.trans_r);
    return *this;
}
//
//  Inversion in self-contained and other form
//
HepTransformation& HepTransformation::invert() {
    trans_r.invert();
    trans_t = trans_r * trans_t.inverse().trans_vec();
    return *this;
}

HepTransformation HepTransformation::inverse() {
    return HepTransformation(trans_r.inverse(), trans_r.inverse() * trans_t.inverse().trans_vec());
}

HepRotation HepTransformation::rot_mat() const {
    return trans_r;
}

HepTranslation HepTransformation::tra_vec() const {
    return trans_t;
}

HepPoint HepTransformation::origin() const {
    Hep3Vector t = trans_t.trans_vec();
    return HepPoint(t.x(), t.y(), t.z());
}

Hep3Vector HepTransformation::unit(int axis) const {
    double i(0.), j(0.), k(0.);
    switch (axis) {
        case 0:
            i = 1.;
            break;
        case 1:
            j = 1.;
            break;
        case 2:
            k = 1.;
            break;
        default:
            k = 1.;
    }
    Hep3Vector n(i, j, k);
    return trans_r * n;
}
// Returns the coordinate of a point in the first system,
//  given its coordinates in the second system
HepPoint HepTransformation::transFrom(const HepPoint& p2) const {
    return trans_r * p2 + trans_t.trans_vec();
}
// Returns the coordinate of a point in the second system,
//  given its coordinates in the first system
HepPoint HepTransformation::transTo(const HepPoint& p1) const {
    return trans_r.inverse() * (p1 - trans_t.trans_vec());
}
// Returns the coordinate of a vector in the first system,
//  given its coordinates in the second system
Hep3Vector HepTransformation::transFrom(const Hep3Vector& v2) const {
    return trans_r * v2;
}
// Returns the coordinate of a vector in the second system,
//  given its coordinates in the first system
Hep3Vector HepTransformation::transTo(const Hep3Vector& v1) const {
    return trans_r.inverse() * v1;
}

Hep3Vector& HepTransformation::transform(Hep3Vector& v) const {
    return v = rot_mat() * v;
}
