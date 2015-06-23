//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchVolElem.hh 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchVolElem
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//	
//
// Copyright Information:
//	Copyright (C) 1997	R. Stroili
//
//------------------------------------------------------------------------

#ifndef DCHVOLELEM_HH
#define DCHVOLELEM_HH

#include <vector>

//----------------------
// Base Class Headers --
//----------------------
#include "DetectorModel/DetVolumeElem.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class DchVolType;
class DchDetector;
class HepPoint;
class Trajectory;
class Transformation;

//		---------------------
// 		-- Class Interface --
//		---------------------

class DchVolElem : public DetVolumeElem {

  //--------------------
  // Instance Members --
  //--------------------

public:

  // Constructors
  DchVolElem(DchVolType* itsType, const char* name, int id,
      const HepTransformation& theAlignment);
  // Copy Constructor
  //    DchVolElem( const DchVolElem& );

  // Destructor
  virtual
  ~DchVolElem();

  // Operators

  // Selectors (const)
  virtual int
  intersect(const Trajectory*, DetIntersection&) const;
  virtual void
  sideIntersect(const Trajectory*, std::vector<DetVolSideIntersection>&,
      double flightDistance, double* arrayOfRange) const;

  bool
  debug(void) const
  {
    return _debug;
  }

  // Modifiers
  //   void createCache();
  void
  setDebug(bool deb)
  {
    _debug = deb;
  }

protected:

  // Helper functions
  void
  setStep(void) const;
  void setStartStepSize(double stepsize);

private:

  bool _extSettedStrtStp;
  double _strtSTEP;

  //  input point is in the BaBar reference frame
  bool
  insideVolume(const HepPoint& point) const;

  // Friends
  friend class DchDetector;

  bool _debug;

};

#endif
