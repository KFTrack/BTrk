//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkParams.cc,v 1.16 2004/08/06 06:31:43 bartoldu Exp $
//
// Description:
//     Implementation of TrkParams class
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Dave Brown, Steve Schaffner
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include <assert.h>
#include "TrkBase/TrkParams.hh"
#include <iostream>
using std::ostream;

//  Construct with the correct dimension
//--------------------------------------
TrkParams::TrkParams(int ndim) :   DifIndepPar(ndim), _weightInverted(false)
//--------------------------------------
{}

//  Construct from an input parameter vector and covariance matrix
//--------------------------------------
TrkParams::TrkParams(const HepVector& pvec,const HepSymMatrix& pcov) :
  DifIndepPar(pvec, pcov), _weightInverted(false)
//--------------------------------------
{}

// copy constructor
//--------------------------------------
TrkParams::TrkParams(const TrkParams& old) :
  DifIndepPar(old), _weightInverted(old._weightInverted), _weightMatrix(old._weightMatrix)
//--------------------------------------
{}

// base implementation of clone function
//--------------------------------------
TrkParams*
TrkParams::clone() const {
//--------------------------------------
  return new TrkParams(*this);
}

//--------------------------------------
TrkParams::~TrkParams() {}
//--------------------------------------

//   Assignment operator
//--------------------------------------
TrkParams& 
TrkParams::operator = (const TrkParams& other){
//--------------------------------------
  if (this == &other) return *this;
  parameter() = other.parameter();
  covariance() = other.covariance();
  _weightInverted =  other._weightInverted ;
  _weightMatrix = other._weightMatrix ;
  return *this;
}

//--------------------------------------
void 
TrkParams::printAll(ostream& os) const {
//--------------------------------------
  DifIndepPar::printAll(os);
}


const HepSymMatrix& TrkParams::weightMatrix() const {

  if (!_weightInverted){
    _weightMatrix = covariance();
    int invStatus;
    _weightMatrix.invert(invStatus);
    _weightInverted = true;
  }
  return _weightMatrix;

}
