#ifndef _HEPPOINT_H_
#define _HEPPOINT_H_

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Vector/ThreeVector.h"

using Point3D = HepGeom::Point3D<double>;

class HepPoint : public Point3D {
   public:
    using Point3D::Point3D;

    HepPoint& operator+=(const CLHEP::Hep3Vector& p) {
        v_[0] += p.x();
        v_[1] += p.y();
        v_[2] += p.z();
        return *this;
    };

    HepPoint& operator-=(const CLHEP::Hep3Vector& p) {
        v_[0] -= p.x();
        v_[1] -= p.y();
        v_[2] -= p.z();
        return *this;
    };
};

#endif
