//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHelixData_001.hh,v 1.6 2007/09/24 21:56:27 gapon Exp $
//
//  Description:
//  Class TrkHelixData_001; a compact representation of a helix and related
//  information.  It is fixed length and contains no pointer data, so it
//  can be included in any persistence scheme.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Infomation;
//	Copyright (C) 2000	Lawrence Berkeley Laboratory
//
// Author(s): Dave Brown 10/24/00
//
//------------------------------------------------------------------------

#ifndef TRKHELIXDATA_HH
#define TRKHELIXDATA_HH

#include "BaBar/BaBarODMGTypes.h"
#include "PDT/PdtPid.hh"
#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/HelixTraj.hh"

template <class T> class ComPackBase;
class ComPackExpFloat;
class ComPackSignedExpFloat;
class ComPackFlatFloat;
class ComPackInt;
class TrkFitSummary;
class BField;
#include "CLHEP/Vector/ThreeVector.h"

class TrkHelixData_001 {
public:
// default constructor for Objectivity
  TrkHelixData_001();
// construct from a HelixTraj 
  TrkHelixData_001(const HelixTraj*);
// construct from packed data
  TrkHelixData_001(const d_ULong params[5], const d_ULong corr[3], const d_ULong& fltlen);
// copy and equivalence are OK
  TrkHelixData_001(const TrkHelixData_001&);
  TrkHelixData_001& operator =(const TrkHelixData_001&);
// NO VIRTUAL DESTRUCTOR.  Yes, this is on purpose, to cut down on overhead
  ~TrkHelixData_001();
// allow separate unpacking of all parameters;
// the helix function _RETURNS OWNERSHIP_
  HelixTraj* helix() const;
  void flightRange(double range[2]) const;
// access to raw data
  const d_ULong& parameters(HelixTraj::ParIndex index) const {
    return _params[index]; }
  const d_ULong& correlations(unsigned icor) const {
    return _corr[icor]; }
  const d_ULong& flightLength() const {
    return _fltlen; }
  void setData(d_ULong params[5],d_ULong corr[3],d_ULong& fltlen);
  int charge(const BField* field) const;
  Hep3Vector momentum(const BField* field,double fltlen) const;
private:
// array of parameter and diagonal covariance terms.
// these are indexed by the same enum defined in HelixTraj.
  d_ULong _params[HelixTraj::NHLXPRM];
// array of data used to store the correlation matrix
  d_ULong _corr[3];
// a single word to store the flightlength range
  d_ULong _fltlen;
// utility functions
  static const ComPackBase<double>& paramPacker(HelixTraj::ParIndex);
  static const ComPackBase<double>& errorPacker(HelixTraj::ParIndex);
  static const ComPackBase<double>& corrPacker(HelixTraj::ParIndex,HelixTraj::ParIndex);
  static void corrTerm(unsigned iterm,
		       HelixTraj::ParIndex& ipar,
		       HelixTraj::ParIndex& jpar);
// unpack the correlation matrix.  The input _must_ be a diagonal matrix,
// hence this function is private
  void correlation(HepSymMatrix& corr,unsigned* ) const;
// statics used for packing
  static const ComPackSignedExpFloat _packd0;
  static const ComPackSignedExpFloat _packz0;
  static const ComPackFlatFloat _packphi0;
  static const ComPackSignedExpFloat _packomega;
  static const ComPackFlatFloat _packlambda;

  static const ComPackExpFloat _packd0err;
  static const ComPackExpFloat _packz0err;
  static const ComPackExpFloat _packphi0err;
  static const ComPackExpFloat _packomegaerr;
  static const ComPackExpFloat _packlambdaerr;

  static const ComPackExpFloat _packcorrfine;
  static const ComPackExpFloat _packcorrcoarse;
  static const ComPackFlatFloat _packcorrflat;
  static const ComPackInt _packsigns;
  static const ComPackSignedExpFloat _packfltlen;
  static const ComPackFlatFloat _packfltrange;
  
  static const double _rescale;
  static const double _mindet;

// unpack a helix with just parameters (no errors!)
  void unpackHelix(HelixTraj&) const;

// persistence
  friend class HelixTrajK;
};

#endif
