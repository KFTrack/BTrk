#ifndef DCHCELL_HH
#define DCHCELL_HH
//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchCell.hh 91 2010-01-14 12:37:23Z stroili $
//
// Description:
//	Class DchCell.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//      M. Kelsey               980129  Add cellID<->wire,layer functions
//
// Copyright Information:
//	Copyright (C) 1997	INFN-Pd
//
//	20020411  Bug fix -- operator<< ought to take const DchCell&
//	20020419  Add protected helixPath(TrkExchangePar...), called from
//		  both TrkFit* and HelixTraj& public interfaces.
//------------------------------------------------------------------------

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <vector>

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "DchGeomBase/DchCellAddr.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class DchFWire;
class DchSWire;
class HepPoint;
class Trajectory;
class TrkFit;
class HelixTraj;
class TrkExchangePar;

class DchCell {

public:

  //  Constructors -- outline only, no sense wire
  DchCell( std::vector<DchFWire*>& );

  //  Copy Constructor
  DchCell( const DchCell& );

  //  constructor from phi rotation of the base cell
  DchCell( const DchCell&, double phi, int cellnum=-1, 
	   const DchSWire* sense=0);

  //  Destructor
  virtual ~DchCell( void ) ;

  //  Operators
  int operator==( const DchCell& ) const;
  //  return, from base cell, a new one "rotating" it by phi 
  //  around the z BaBar axis
  //   DchCell& rotate(double phi);

  //  returns path length of track in cell, it needs a starting flt
  float intersect( const Trajectory*, double range[2] ) const;
  float path( const Trajectory*, double trkRange[2], double wireFlt=0 ) const;

  //  returns path length of track in cell and intersection points, it needs 
  //  a starting flt
  float intersect( const Trajectory*, double range[2], 
		   std::vector<HepPoint>& vec ) const;
  float path( const Trajectory*, double trkRange[2], double wireFlt, 
              std::vector<HepPoint>& vec ) const;

  //  fast path-length algorithm which uses local helix parameters
  //  NOTE:  See below for internal interface
  double helixPath(const TrkFit* traj, double fltlen) const;
  double helixPath(const HelixTraj& theHelix, double fltlen) const;

  //   NOT implemented. Will ever be used? 
//   HepPoint intersectUp( const Trajectory*, double range[2] ) const;
//   HepPoint intersectDown( const Trajectory*, double range[2] ) const;

  void print(std::ostream& o=std::cout) const ;
  void printAll(std::ostream& o=std::cout) const ;

  //  Selectors (const)
  int cellNum(void) const  { return _cellNum; }
  int wireNum(void) const  { return DchCellAddr::wireIs(_cellNum); }
  int layerNum(void) const { return DchCellAddr::layerIs(_cellNum); }

  const DchSWire* getSense(void) const { return _sense; }

  bool inCell(const HepPoint& p ) const ;

  const std::vector<DchFWire*>& cellWires() const { return _cellWires; }

protected:
  double helixPath(const TrkExchangePar& helix, const HepPoint& ref,
		   const HepPoint& hit, double fltlen) const;

private:

  //  Data members
  std::vector<DchFWire*> _cellWires;		// field wires delimiting the 
						// drift cell
  const DchSWire* _sense;			// sense wire for cell
  int _cellNum;					// cell number
};

std::ostream& operator << (std::ostream& o, const DchCell&);

#endif // DCHCELL_HH
