//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchGasSet.hh 91 2010-01-14 12:37:23Z stroili $
//
// Description:
//	Class DchGasSet
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	<R. Stroili>		<originator>
//	
//
// Copyright Information:
//	Copyright (C) 1994	<INFN & Universita' di Padova>
//
//------------------------------------------------------------------------

#ifndef DCHGASSET_HH
#define DCHGASSET_HH

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------
#include <vector>

//----------------------
// Base Class Headers --
//----------------------
#include "DetectorModel/DetSet.hh"

//		---------------------
// 		-- Class Interface --
//		---------------------

//template<class X> 
class DchGasSet : public DetSet {

//--------------------
// Instance Members --
//--------------------

public:

  //  Constructors
  DchGasSet(char* name, int iset) : DetSet(name,iset) {;}

  //  Destructor
  virtual ~DchGasSet( );

  //  Operators
//
//  Find the next element intersected by a given Trajectory from a given 
//  starting intersection
//  Don't use this to loop over elements, as it is inefficient.
//
  bool nextIntersection( const Trajectory* traj,DetIntersection& next,
			 DetIntersection* previous = 0 ) const;
//
//  Make a complete list of interesected DetElems for a given trajectory.  
//  The (optional) supplied range overrides the trajectory range
//
   void intersection( std::vector<DetIntersection>&,
 		     const Trajectory*,double* myrange=0,
 		     bool clear=true ) const;

private:
  //  Copy Constructor
  DchGasSet( const DchGasSet& );

};

#endif
