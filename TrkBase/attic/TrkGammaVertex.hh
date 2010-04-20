//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkGammaVertex.hh,v 1.10 2007/08/29 19:13:58 ttanabe Exp $
//
// Description:
//      A vertexing class specialized in finding photon conversions.
//      The algorithm combines POCA and parallelism in a least-squares fit.
//      Intended for use in the photon conversion list at the Beta-level
//      (SimpleComposition) and at the track-level as a V0 finder (TrkFixup).
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Tomohiko Tanabe, David Brown, LBNL (12/14/2006)
//
// Copyright Information:
//      Copyright (C) 2006              Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------
#ifndef TrkGammaVertex_HH__
#define TrkGammaVertex_HH__

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkVertex.hh"
#include "TrkBase/TrkErrCode.hh"

class TrkRecoTrk;
class TrkFit;
class TrkDifTraj;
class BbrPointErr;
class PdtEntry;
class HelixTraj;
class BField;
class HepHistogram;
class HepTuple;
class HepTupleManager;

// Class interface //
class TrkGammaVertex : public TrkVertex {
public:
  // constructors
  TrkGammaVertex();
  ~TrkGammaVertex();
  enum AlgoMode { Conversion, Dalitz };

  /* This class takes two tracks as input, vertex them using
     the algorithm assuming they come from a photon conversion.
     This is described in some detail in the .cc file.

     Depending on the information you already have about the two tracks,
     you can invoke three different constructors which are described below.
   */

  /* Use this constructor if you don't know much about the two tracks.
     This is probably the constructor you want to use by default. */
  TrkGammaVertex(const TrkRecoTrk* trk1, const TrkRecoTrk* trk2,
                 enum AlgoMode algomode);
  
  /* Use this constructor if the track pair is already vertexed
     (by vertexers like TrkPocaVertex) */
  TrkGammaVertex(const TrkVertex* vtx, enum AlgoMode algomode);

  /* copy & clone constructors, etc */
  TrkGammaVertex(const TrkGammaVertex& other);
  TrkGammaVertex& operator = (const TrkGammaVertex& other);
  virtual TrkGammaVertex* clone() const;

  // accessors
  double flt1() const { return _flt1; }
  double flt2() const { return _flt2; }
  double doca() const { return _doca; }
  const HepMatrix& xpCov() const { return _vxp; }
  const HepHistogram* fillChisq(
    HepTupleManager* manager, const char* name,
    HepTuple* tuple, const HepPoint& tpos,
    int gind1, int gind2) const;

  const HepVector& gpar() const { return _gpar; }
  const HepSymMatrix& gcov() const { return _gcov; }

private:
  // data members
  double _flt1;
  double _flt2;
  double _doca;
  HepPoint _pp; // parallel point
  AlgoMode _algomode;
  HepMatrix _vxp; // position-momentum correlation matrix
  HelixTraj* _htrajs[2];
  HepVector _gpar;
  HepSymMatrix _gcov;

// indices for gamma parameters
  enum gind{om1=0,om2,phi0,dip,d0,z0,dphi0,ddip,tf,ngind};
  enum rind{xind=0,yind,zind};

  // private methods: algorithm implementations and utility methods
  void vertexIt(
      const TrkRecoTrk* const & trk1,
      const TrkRecoTrk* const & trk2);

  TrkErrCode transParallel(const TrkDifTraj** intrajs, double* flt, HepPoint& pp) const;
  
  bool fitSelect(const TrkFit* fit1,const TrkFit* fit2) const;

  void fillParameters(const HepVector& gpar, HepVector* trkpars) const;

  TrkErrCode minimize(HelixTraj** htrajs,HepVector& gpar,HepSymMatrix& gcov,double& chisq) const;
  void findVertex(const BField& bfield, const HepVector& gpar,
      const HepSymMatrix& gcov, const HepPoint& pp);

  // these are used for caching PdtEntry lookups
  static const PdtEntry * getGammaPdt();
  static const PdtEntry * getElectronPdt();
  // constants used by the algorithm
  static const double _maxdoca;
  static const double _maxddip;
  static const double _awconv;
  static const double _awdalitz;
};

#endif

