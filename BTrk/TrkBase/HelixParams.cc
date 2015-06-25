//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: HelixParams.cc,v 1.19 2004/08/06 06:31:40 bartoldu Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Steve Schaffner
//------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BbrGeom/BbrAngle.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include <assert.h>
#include <iostream>
using std::endl;
using std::ostream;
using namespace CLHEP;
//----------------------------------------------------------------------
HelixParams::HelixParams(const HepVector& params, const HepSymMatrix& pcov) 
  :  TrkParams(params,pcov) { 
//----------------------------------------------------------------------
    assert(parvec.num_row() == NHLXPRM);
    assert(parcov.num_row() == NHLXPRM);
    parvec[phi0Index] = BbrAngle(parvec[phi0Index]).rad();
}

//----------------------------------------------------------------------
HelixParams::~HelixParams(){}

//----------------------------------------------------------------------
// Output functions: simple, detailed, and stream

void HelixParams::print(ostream& o) const {
  o << "d0="  << d0() << " phi0="   << phi0() << " omega=" << omega()
    << " z0=" << z0() << " tanDip=" << tanDip();
}

void HelixParams::printAll(ostream& o) const {
  print(o);
  o << endl << covariance() << endl;
}

ostream& operator<<(ostream& o, const HelixParams& helix) {
  helix.print(o);
  return o;
}

