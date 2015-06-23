//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id:  $
//
// Description:
//	Class DchPhiSegmVolElem
//
// Environment:
//
//
// Author List:
//	
//
//
// Copyright Information:
//
//
//------------------------------------------------------------------------

#ifndef DCHPHISEGMVOLELEM_HH
#define DCHPHISEGMVOLELEM_HH

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

class DchPhiSegmVolElem : public DetVolumeElem {

  //--------------------
  // Instance Members --
  //--------------------

public:

  // Constructors
  DchPhiSegmVolElem(DchVolType* itsType, const char* name, int id,
      const HepTransformation& theAlignment);
  // Copy Constructor
  //    DchPhiSegmVolElem( const DchPhiSegmVolElem& );

  // Destructor
  virtual
  ~DchPhiSegmVolElem();

  // Operators

  // Selectors (const)
  virtual int
  intersect(const Trajectory*, DetIntersection&) const;

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

private:

  //  input point is in the BaBar reference frame
  bool
  insideVolume(const HepPoint& point, double &extraDelta, int nRot=-1) const;

  bool
  chckIntrsctInAglLmts(const Trajectory*, DetIntersection&) const;

  // Friends
  friend class DchDetector;

  bool _debug;

};

#endif
