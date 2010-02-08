//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHelixData.hh,v 1.14 2003/08/07 00:17:04 brownd Exp $
//
//  Description:
//  Class TrkHelixData; a compact representation of a helix and related
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
#include "TrkBase/TrkEnums.hh"

template <class T> class ComPackBase;
class ComPackExpFloat;
class ComPackSignedExpFloat;
class ComPackFlatFloat;
class ComPackInt;
class TrkFitSummary;
class BField;

class TrkHelixData {
public:
// default constructor for Objectivity
  TrkHelixData();
// construct from a HelixTraj plus chisq probability and PID information
  TrkHelixData(const HelixTraj*,
	       double chiprob,PdtPid::PidType hypo,
	       TrkEnums::PackFlag flag=TrkEnums::KalFit);
// construct from packed data
  TrkHelixData(d_ULong params[5], d_ULong corr[3], d_ULong fltlen);
// copy and equivalence are OK
  TrkHelixData(const TrkHelixData&);
  TrkHelixData& operator =(const TrkHelixData&);
// NO VIRTUAL DESTRUCTOR.  Yes, this is on purpose, to cut down on overhead
  ~TrkHelixData();
// allow separate unpacking of all parameters;
// the helix function _RETURNS OWNERSHIP_
  HelixTraj* helix() const;
  double chisqProb() const;
  void flightRange(double range[2]) const;
  PdtPid::PidType pidType() const;
  TrkEnums::PackFlag fitFlag() const;
// return a single object
  TrkFitSummary fitSummary() const;
// access to raw data
  const d_ULong& parameters(HelixTraj::ParIndex index) const {
    return _params[index]; }
  const d_ULong& correlations(unsigned icor) const {
    return _corr[icor]; }
  const d_ULong& flightLength() const {
    return _fltlen; }
  int charge(const BField* field) const;
private:

  friend class KalMiniTrkK;
  friend class TrkHelixDataK;
  friend class TrkKalMiniCompositeK;

// array of parameter and diagonal covariance terms.
// these are indexed by the same enum defined in HelixTraj.
  d_ULong _params[HelixTraj::NHLXPRM];
// array of data used to store the correlation matrix, plus
// a few odd bits to store PID and a usage flag.
  d_ULong _corr[3];
// a single workd to store the flightlength range and chisquared probability
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
  static const ComPackInt _packpid;
  
  static const ComPackExpFloat _packchiprob;
  static const ComPackSignedExpFloat _packfltlen;
  static const ComPackFlatFloat _packfltrange;
  static const ComPackInt _packflag;

  static const double _rescale;
  static const double _mindet;
};

#endif
