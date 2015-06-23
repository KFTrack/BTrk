//
//specific to this code, taken from framework

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/HepPoint.h"
#include "mu2eFast/mu2eDSField.hh"
#include "Framework/AppFileName.hh"
#include "PacEnv/PacConfig.hh"

#include <iostream>

int main(int argc, char* argv[]) {
  
  gconfig.verbose(true);
  gconfig.parseargs(argc, argv);
  double dfactor = gconfig.getfloat("distortionfactor",1.0);
  
  mu2eDSField dsfield("mu2eFast/DSMap_unfmt_rad100.dat",dfactor);
  
  double bnom = dsfield.bFieldNominal();
  std::cout <<" nominal field = " << bnom << std::endl;

  // 
  // just something to scan out the z-axis +- 5 meters; steps are mm.  I've
  // checked a few individual test points by hand.

  std::cout << "z scan " << std::endl;
  double x = 0.;
  double y = 0.;
  for (int i = -200; i <=200; i+=10){
    double z = i;
    HepPoint testpoint = HepPoint(x,y,z);
    Hep3Vector bfield = dsfield.bFieldVect(testpoint);
    std::cout << "testpoint is " << testpoint << " and BField is " << bfield << std::endl;
  }
  
  std::cout << "micro z scan " << std::endl;
  x = 50.;
  y = 0.;
  for (int i = -50; i <=50; i+=1){
    double z = i;
    HepPoint testpoint = HepPoint(x,y,z/10.0);
    Hep3Vector bfield = dsfield.bFieldVect(testpoint);
    std::cout << "testpoint is " << testpoint << " and BField is " << bfield << std::endl;
  }
  
  
  std::cout << "r scan" << std::endl;
  double z=0.0;
  for(int i=0;i<90;i++){
    x = i;
    HepPoint testpoint = HepPoint(x,y,z);    
    Hep3Vector bfield = dsfield.bFieldVect(testpoint);
    std::cout << "testpoint is " << testpoint << " and BField is " << bfield << std::endl;
  }

  std::cout << "phi scan" << std::endl;
  double rho=90.0;
  for(int i=0;i<100;i++){
    double phi = i*6.283185/100;
    x = rho*cos(phi);
    y = rho*sin(phi);
    HepPoint testpoint = HepPoint(x,y,z);    
    Hep3Vector bfield = dsfield.bFieldVect(testpoint);
    std::cout << "testpoint is " << testpoint << " and BField is " << bfield << std::endl;
  }
  
  return 0;
}