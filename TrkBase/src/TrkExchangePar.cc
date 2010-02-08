//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkExchangePar.cc,v 1.19 2004/08/06 06:31:40 bartoldu Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Steve Schaffner
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "BbrGeom/BbrAngle.hh"
#include <assert.h>
#include <iostream>
using std::endl;
using std::ostream;

//----------------------------------------------------------------------
TrkExchangePar::TrkExchangePar(const HepVector& inV, const HepMatrix& inErr) 
  : paramVec(inV) { 
//----------------------------------------------------------------------
    assert(paramVec.num_row() == nParam);
    assert(inErr.num_row() == nParam);
    assert(inErr.num_col() == nParam);
    paramErr.assign(inErr);
    paramVec[ex_phi0] = BbrAngle(paramVec[ex_phi0]).rad();
}

//----------------------------------------------------------------------
TrkExchangePar::TrkExchangePar(const HepVector& inV, const HepSymMatrix& err) 
  :  paramVec(inV), paramErr(err) { 
//----------------------------------------------------------------------
    assert(paramVec.num_row() == nParam);
    assert(paramErr.num_row() == nParam);
    paramVec[ex_phi0] = BbrAngle(paramVec[ex_phi0]).rad();
}

//----------------------------------------------------------------------
TrkExchangePar::TrkExchangePar(const HepVector& inV) 
  : paramVec(inV), paramErr(nParam, 1) {
//----------------------------------------------------------------------
    assert(paramVec.num_row() == nParam);
    paramVec[ex_phi0] = BbrAngle(paramVec[ex_phi0]).rad();
}

//----------------------------------------------------------------------
TrkExchangePar::TrkExchangePar(double d0In, double phi0In, double omegaIn, 
			       double z0In, double tanDipIn) 
  : paramVec(nParam), paramErr(nParam, 1) { 
//----------------------------------------------------------------------
    paramVec[ex_d0]     = d0In;
    paramVec[ex_phi0]   = phi0In;
    paramVec[ex_omega]  = omegaIn;
    paramVec[ex_z0]     = z0In;
    paramVec[ex_tanDip] = tanDipIn;
    paramVec[ex_phi0] = BbrAngle(paramVec[ex_phi0]).rad();
}

//----------------------------------------------------------------------
TrkExchangePar::~TrkExchangePar(){}

//----------------------------------------------------------------------
// Output functions: simple, detailed, and stream

void TrkExchangePar::print(ostream& o) const {
  o << "d0="  << d0() << " phi0="   << phi0() << " omega=" << omega()
    << " z0=" << z0() << " tanDip=" << tanDip();
}

void TrkExchangePar::printAll(ostream& o) const {
  print(o);
  o << endl << covariance() << endl;
}



ostream& operator<<(ostream& o, const TrkExchangePar& helix) {
  helix.print(o);
  return o;
}

