//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkRep.hh,v 1.69 2007/09/24 21:56:27 gapon Exp $
//
// Description: Base class for internal track representation classes -- e.g. 
//   HelixRep, KalRep.  Owns and maintains a TrkHitVector; and keeps a 
//   pointer to the track that owns the Rep.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Steve Schaffner
//
//------------------------------------------------------------------------

#ifndef TRKREP_HH
#define TRKREP_HH
// BTrk includes
#include "BTrk/TrkBase/TrkDirection.hh"
#include "BTrk/TrkBase/TrkFitStatus.hh"
#include "BTrk/TrkBase/TrkFit.hh"
#include "BTrk/TrkBase/TrkHitUpdater.hh"
#include "BTrk/TrkBase/TrkHit.hh"
#include "BTrk/TrkBase/TrkT0.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
// general includes
#include "CLHEP/Vector/ThreeVector.h"
#include <memory>
#include <vector>
#include <functional>
#include <iosfwd>
#include <memory>

class TrkErrCode;
class HepPoint;
class ChisqConsistency;
class BbrPointErr;
class BbrVectorErr;

//namespace BTrk {
// A TrkHitVector is just a vector of pointers to hits.
typedef std::vector<TrkHit*> TrkHitVector;
// hits are sorted by flightlength
struct hitsort : public std::binary_function<TrkHit*, TrkHit*, bool> {
  bool operator()(TrkHit* x, TrkHit* y) { 
    return  x->fltLen() < y->fltLen();
  }
};

// Class interface //
class TrkRep : public TrkFitStatus, public TrkFit, public TrkHitUpdater {

  public:
    // construct from a hitlist; the rep takes ownership of the hits
    TrkRep(TrkHitVector const& inTrkHits, TrkParticle const& tpart );
    virtual ~TrkRep();
    TrkRep&   operator = (const TrkRep&) = delete;
    bool operator== (const TrkRep&);

    //******************************************
    // Global quantities:
    //******************************************
    virtual ChisqConsistency    chisqConsistency() const;
    virtual int               nActive()      const;
    double                    startValidRange() const;
    double                    endValidRange()   const;
    virtual double            startFoundRange() const;
    virtual double            endFoundRange()   const;
    //  t0
    TrkT0 const& t0() const { return _trkt0; }
    // flightlength associated with t0.
    double flt0() const { return _flt0; }
    TrkParticle const& particleType() const { return _tpart; }
    // for now, allow setting t0.  Must be coherent with flt0
    void setT0(TrkT0 const& t0,double flt0) { _trkt0 = t0; _flt0 = flt0;}

    //******************************************
    // Information about track at a given position
    //******************************************
    virtual HepPoint                  position(double fltL)         const;
    virtual CLHEP::Hep3Vector                direction(double fltL)        const;
    virtual double                    arrivalTime(double fltL)      const;
    virtual BbrPointErr               positionErr(double fltL)      const;
    virtual BbrVectorErr              directionErr(double fltL)     const;

    //******************************************
    // Hit (list) handling
    //******************************************
    // Simple implementations of these are present in the base class; 
    //   complicated reps (e.g. Kalman) may wish to override.
    virtual void		    addHit(TrkHit *theTrkHit);
    virtual void		    removeHit(TrkHit *theTrkHit);
    virtual void		    activateHit(TrkHit *theTrkHit);
    virtual void		    deactivateHit(TrkHit *theTrkHit);
    virtual TrkHitVector const&     hitVector() const {return _hitvec;}
    virtual TrkHitVector&	    hitVector() {return _hitvec;}
    virtual void		    updateTrkHits();
    virtual bool		    resid(const TrkHit *theTrkHit,
	double &residual, double &residErr,
	bool exclude=false) const;

    //****************************************** 
    // Fitting stuff
    //******************************************
    virtual TrkErrCode fit() = 0;

  protected:
    // new data members for mu2e
    TrkT0 _trkt0; // t0 value
    double _flt0; // flightlength associated with t=t0.
    TrkParticle _tpart; // particle type of this rep
    TrkHitVector  _hitvec; // hits used in this rep
    void sortHits(); // allow subclasses to sort hits
};

#endif
