#ifndef _HEPPOINT_H_
#define _HEPPOINT_H_

#include <iostream>
#include "CLHEP/Geometry/BasicVector3D.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Vector/ThreeVector.h"

using Point3D = HepGeom::Point3D<double>;

class HepPoint : public Point3D {
   public:
    HepPoint() : Point3D() { std::cout << "CTor 1" << std::endl; }

    HepPoint(double x1, double y1, double z1) : Point3D(x1, y1, z1) {
        std::cout << "CTor 2" << std::endl;
    }

    HepPoint(const double* a) : Point3D(a) { std::cout << "CTor 3" << std::endl; }

    HepPoint(const Point3D<double>& v) : Point3D(v) { std::cout << "CTor 4" << std::endl; }

    HepPoint(const CLHEP::Hep3Vector& v) : Point3D(v) { std::cout << "CTor 5" << std::endl; }

    HepPoint(const HepGeom::BasicVector3D<double> v) : Point3D(v) {
        std::cout << "CTor 6" << std::endl;
    }

    HepPoint(const HepGeom::BasicVector3D<float> v) : Point3D(v) {
        std::cout << "CTor 7" << std::endl;
    }

    HepPoint& operator+=(const CLHEP::Hep3Vector& p) {
        v_[0] += p.x();
        v_[1] += p.y();
        v_[2] += p.z();
        return *this;
    }

    HepPoint& operator-=(const CLHEP::Hep3Vector& p) {
        v_[0] -= p.x();
        v_[1] -= p.y();
        v_[2] -= p.z();
        return *this;
    }
};

#endif
