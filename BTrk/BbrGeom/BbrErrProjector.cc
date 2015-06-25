//$Id: BbrErrProjector.cc 491 2010-01-13 16:59:16Z stroili $
//Implementation of BbrErrProjector
//A. Snyder
//SLAC  Fri Jul 16 14:51:01 PDT 1999
#include "BTrk/BaBar/BaBar.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "BTrk/BbrGeom/BbrErrProjector.hh"
using namespace CLHEP;

BbrErrProjector::Code
BbrErrProjector::localError(double sigma,
			    const Hep3Vector& direction,
			    const Hep3Vector& uHat,
			    const Hep3Vector& vHat,
			    HepSymMatrix& error) {

  if(error.num_row()!=2) return matrixWrongSize;

  double zeta=uHat*vHat;
  if(fabs(zeta)<0.001) {
    Hep3Vector nHat=uHat.cross(vHat);
    double an=direction*nHat;
    if(fabs(an)<0.00001) return trackParallel;
    double au=direction*uHat;
    double av=direction*vHat;
    error(1,1)=au*au+an*an;
    error(1,2)=au*av;
    error(2,2)=av*av+an*an;
    error*=sigma*sigma/(an*an);
    return ok;

  }else if(fabs(zeta)<0.99) {
    Hep3Vector uHatP=(uHat-vHat*zeta);
    uHatP*=1/sqrt(1-zeta*zeta);
    HepSymMatrix ep(2);
    Code ecode=localError(sigma,direction,uHatP,vHat,ep);
    if(ecode!=ok) return ecode;
    error(1,1)=ep(1,1)*(1-zeta*zeta);
    error(1,1)+=ep(2,2)*zeta*zeta;
    error(1,1)+=2*ep(1,2)*sqrt(1-zeta*zeta)*zeta;
    error(2,2)=ep(2,2);
    error(1,2)=ep(2,2)*zeta;
    error(1,2)+=ep(1,2)*sqrt(1-zeta*zeta);
    return ok;
  }else {
    return uvTooClose;
  }

}
