// -*- C  ++ -*-
// $Id: Translation.cc 192 2009-03-04 12:20:53Z stroili $
// ---------------------------------------------------------------------------
//
// This file is part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the HepTranslation class. 
//
// .SS History
// Author Victoria Novotny (LBL)

#ifdef __GNUC__
#pragma interface
#endif

#include "CLHEP/Geometry/Translation.h"

HepTranslation::~HepTranslation() {}

HepTranslation::HepTranslation(): vec(0.0, 0.0, 0.0) {}

HepTranslation::HepTranslation(const Hep3Vector & m):
    vec(m) {}

HepTranslation::HepTranslation(const HepTranslation & m){ 
    vec = m.trans_vec();
}


HepTranslation & HepTranslation::operator = (const HepTranslation &m) {
    vec = m.trans_vec();
    return *this;
}

double HepTranslation::x() const { return vec.x(); }
double HepTranslation::y() const { return vec.y(); }
double HepTranslation::z() const { return vec.z(); }

Hep3Vector HepTranslation::trans_vec() const { return vec; }

HepTranslation HepTranslation::inverse() const {
    HepTranslation m(-vec);
    return m; 
}

HepTranslation & HepTranslation::invert() {
    return *this=inverse();
}

HepTranslation & HepTranslation::translateX(double p) {
    vec.setX(vec.x()+p);
    return *this;
}

HepTranslation & HepTranslation::translateY(double p) {
    vec.setY(vec.y()+p);
    return *this;
}

HepTranslation & HepTranslation::translateZ(double p) {
    vec.setZ(vec.z()+p);
    return *this;
}

HepTranslation & HepTranslation::operator += (const HepTranslation & p) {
    vec += p.trans_vec();
    return *this;
}

Hep3Vector & HepTranslation::transform(Hep3Vector & v) const {
    return v;
}

HepTranslation operator + (const HepTranslation & q, const HepTranslation & p){
    HepTranslation r(Hep3Vector(q.trans_vec() + p.trans_vec()));
    return r;
}

