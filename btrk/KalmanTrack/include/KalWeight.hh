// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalWeight.hh,v 1.13 2008/03/17 13:14:30 brownd Exp $
//
//  Description: KalWeight
//  define the weight aspect of a KalmanSite: this is just an information vector
//  and an information matrix
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 8/22/97
//------------------------------------------------------------------------------
#ifndef KALWEIGHT_HH
#define KALWEIGHT_HH
#include <iostream>
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "TrkBase/TrkSimpTraj.hh"
class KalParams;

class KalWeight{
public:
  KalWeight();
  KalWeight(int); // specify the dimension
  KalWeight(const HepVector&, const HepSymMatrix&);
  KalWeight(const KalParams&);
  KalWeight(const KalWeight&);
  ~KalWeight(){;}
  unsigned nPar() const { return _wvector.num_row(); }
  const HepVector& weightVector() const { return _wvector; }
  const HepSymMatrix& weightMatrix() const { return _wcov; }
  KalWeight& operator = (const KalWeight& other);
  KalWeight& operator = (const KalParams& other);
// addition operation, used in site processing
  KalWeight& operator += (const KalWeight& other);
// test matrix status after inversion
  bool matrixOK() const { return _status == 0; }
  int status() const { return _status; }
// invalidate to force re-computation
  void invalidate() { _status = 1;}
// printout
  void printAll(std::ostream& os=std::cout) const {
    os << "Weight vector " << _wvector;
    os << "Weight matrix " << _wcov;
  }
  void print(std::ostream& os=std::cout) const {
    os << "Weight vector " << _wvector;
    os << "Weight matrix " << _wcov;
  }
private:
  HepVector _wvector;
  HepSymMatrix _wcov;
  int _status;
};
#endif
