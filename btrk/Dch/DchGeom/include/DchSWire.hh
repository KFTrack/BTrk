//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchSWire.hh 91 2010-01-14 12:37:23Z stroili $
//
// Description:
//	Class DchSWire.
//      Dch sense wire class
//      
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//	
//
// Copyright Information:
//	Copyright (C) 1997	INFN & Padova University
//
//------------------------------------------------------------------------

#ifndef DCHSWIRE_HH
#define DCHSWIRE_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CLHEP/Geometry/HepPoint.h"
#include "DchGeom/DchSagTraj.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class HepTransformation;
class DchLayer;
class DchWirePar;
//class DetAlignElem;
#include <iosfwd>

//		---------------------
// 		-- Class Interface --
//		---------------------

class DchSWire {

//--------------------
// Declarations     --
//--------------------

    // Typedefs, consts, and enums

//--------------------
// Instance Members --
//--------------------

public:

    // Constructors
  DchSWire(const HepPoint& rP, const HepPoint& fP, double sag=0);

    // Copy Constructor
//    DchSWire( const DchSWire& );

    // Destructor
  virtual ~DchSWire( );

    // Operators
    
//    DchSWire&       operator= ( const DchSWire& );
//    virtual int operator==( const DchSWire& ) const;
//            int operator!=( const DchSWire& ) const;

    // Selectors (const)
  double getSag(void) const { return _sag; }
  const HepPoint* getRearPoint(void) const { return &_rear; }
  const HepPoint* getForwPoint(void) const { return &_forward; }
  const DchSagTraj* getTraj(void) const { return &_traj; }
  double xRear(void) const { return _rear.x(); }
  double yRear(void) const { return _rear.y(); }
  double zRear(void) const { return _rear.z(); }
  double xForw(void) const { return _forward.x(); } 
  double yForw(void) const { return _forward.y(); } 
  double zForw(void) const { return _forward.z(); } 
  double xMid(void) const { return rMid() * cos(phi()); }
  double yMid(void) const { return rMid() * sin(phi()); }
  double rEnd(void) const { return _rend; }
  const DchLayer* layer(void) const { return _lay; }
  int Id(void) const { return _id; }
  double zLength(void) const { return getForwPoint()->z()-
				       getRearPoint()->z(); }
  double dPhiz(void) const { return _twist; }
  double zEndDC(void) const { return 0.5*zLength(); }
  double rMid(void) const { return rEnd() * cos( dPhiz() ); }
  double phiE(void) const { return _phiend; }
  double phi(void) const { return _phiend + _twist; }
  double stereo(void) const 
                        { return 2.*rEnd()*sin(dPhiz())/zLength(); }
  double stDip(void) const { return rEnd() - rMid(); }
  // in the local DCH reference
  double radiusDC(double z) const { return rEnd() - stDip() * (1.-  
					   z*z/(zEndDC()*zEndDC())); } 
  double phiDC(double z) const { return phi() + dPhizDC(z); }
  double xWireDC(double z) const { return radiusDC(z)*cos(phiDC(z)); }
  double yWireDC(double z) const { return radiusDC(z)*sin(phiDC(z)); }
  double dPhizDC(double z) const { return atan2( z*stereo(),rMid() ); }

  // direction axis defining the local coordinate system of the sense
  // wire (eventually needed by calibration)
  Hep3Vector yAxis( double z=0. );
  const Hep3Vector& zAxis( void ) const { return _traj.rawDirection(); }

  // alignment member functions
//  void wireAlign( const DchWirePar& corr );   // single wire correction
//  void wireAlign( const DetAlignElem& align );  // DetAlignElem
                                                      // correction
//  void wireAlign( const HepTransformation& transf );  // global transform 
                                                      // correction
//  void wireAlign( const HepTransformation& reartransf,   // plates misalignment
//		  const HepTransformation& forwtransf ); // correction

  void print(std::ostream &o) const;
  void printInfo(std::ostream &o) const;
    // Modifiers

protected:

  void setLayerPtr(const DchLayer* lay) { _lay = lay; }

private:

    // Friends
  friend class DchLayer;

  // Data members: THESE ARE IN THE GLOBAL (BaBar) REFERENCE SYSTEM
  DchSagTraj      _traj;  // wire trajectory
  HepPoint        _rear;  // wire position on rear end-plate
  HepPoint     _forward;  // wire position on forward end-plate
  double           _sag;  // wire sagitta
  const DchLayer*  _lay;  // pointer to layer to which the wire belongs
  double          _rend;  // radius at rear end plate
  double        _phiend;  // phi angle at rear end plate
  double         _twist;  // twist angle between mid and rear chamber
  int               _id;  // sense wire (cell) identifier
};

std::ostream& operator << (std::ostream& o, const DchSWire&);

#endif // DCHSWIRE_HH
