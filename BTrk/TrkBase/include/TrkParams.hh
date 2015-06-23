//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkParams.hh,v 1.19 2004/08/06 06:31:43 bartoldu Exp $
//
// Description:
//  Class to the parameters of a track.  This class should be thought
//  of as a structure, storing the data objects which together form the
//  parameter set of a track. The
//  direct use of this class will be to describe 'simple' trajectories 
//  (5-parameter helix, parabolic tracks, 6 parameter helix, 
//  even 9 parameter helix).  
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
//--------------------------------------------------------------------------
#ifndef TRKPARAMS_HH
#define TRKPARAMS_HH

#include "difAlgebra/DifIndepPar.hh"

#include <iosfwd>

class TrkParams : public DifIndepPar {
public:
  //**************************
  //  Constructors and such
  //**************************
  TrkParams(int nParam);
  //  Construct from a vector and error matrix
  TrkParams(const HepVector&,const HepSymMatrix&);
  //  Copy constructor
  TrkParams(const TrkParams&);
  TrkParams& operator = (const TrkParams&);
  virtual TrkParams* clone() const;
  virtual ~TrkParams();

  //**************************
  //  access
  //**************************
  void printAll(std::ostream& os) const;
  
  // overloaded to invalidate caching. dangerous because non-virtual
  inline HepSymMatrix& covariance() ;
  inline const HepSymMatrix& covariance() const { return DifIndepPar::covariance(); }

  const HepSymMatrix& weightMatrix() const ;

private:
  mutable bool _weightInverted; //is the weight matrix cached?
  mutable HepSymMatrix _weightMatrix; //the inverse of the weight matrix
};

inline HepSymMatrix& TrkParams::covariance()
{
  _weightInverted=false ; 
  return DifIndepPar::covariance();
}

#endif

