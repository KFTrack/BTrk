//
// BaBar includes
#include "BaBar/BaBar.hh"
#include "mu2eFast/mu2eGradientField.hh"

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <limits>
#include <cmath>
#include <assert.h>

using namespace std;

mu2eGradientField::mu2eGradientField(double b0,double z0, double b1, double z1, double rmax) :
_b0(0,0,b0),_b1(0,0,b1),_z0(z0),_z1(z1),_rmax(rmax)
{
  assert(z0 < z1 && z1 < 0.0);
// compute the gradient
  _grad = (b1-b0)/(z1-z0);
}

mu2eGradientField::~mu2eGradientField(){}

// BaBar interface.  Note we have to change units here to the BaBar conventions
Hep3Vector
mu2eGradientField::bFieldVect (const HepPoint &point)const {
  if(point.z() < _z0)
    return _b0;
  else if(point.z() > _z1)
    return _b1;
  else {
    double bgrad = _grad*(point.z()-_z0);
// work in cylindrical coordinates
    double bz = _b0.z()+bgrad;
    double rad = point.perp();
    double bx = -_grad*point.x()/rad;
    double by = -_grad*point.y()/rad;
    return Hep3Vector(bx,by,bz);
  }
}
