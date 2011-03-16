//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchGasSet.cc 91 2010-01-14 12:37:23Z stroili $
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
//	Copyright (C) 1996	<INFN & Universita' di Padova>
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeom/DchGasSet.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CLHEP/Utilities/CLHEP.h"
#include <vector>
using std::cout;
using std::endl;

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------
static const double _epsilon = 1.0e-4;

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
// DchGasSet::DchGasSet(char* name,int iset) 
//   : DetSet(name,iset) 
// {
// }

//--------------
// Destructor --
//--------------
DchGasSet::~DchGasSet() 
{;}

//-------------
// Methods   --
//-------------
    
//-------------
// Operators --
//-------------
//
//  Navigation; find the 'next' element, given a trajectory and a current
//  path length.  This is the stupidest imaginable implementation, trying
//  to force an intersection with every element, and testing on the return
//  value.  Smarter sets should overwrite this function!!!
//
bool 
DchGasSet::nextIntersection( const Trajectory* traj,
                             DetIntersection& nextinter,
			     DetIntersection* lastinter ) const 
{
  cout  << "nextIntersection" <<endl;
  //
  //  Make sure the set is ready
  //
  if ( !isReady() ) makeReady();
  //
  //  Setup the limiting intersections; use the input last intersection if 
  //  provided, else start from the begining of the trajectory range
  //
  bool intersected = false;
  DetIntersection mininter(0,traj,traj->hiRange());
  DetIntersection current(0,traj,traj->lowRange(),traj->hiRange());

  if ( lastinter ) current.pathlen = lastinter->pathrange[1]+_epsilon;
  //
  //  Loop over the contained elements, and intersect the trajectory with them
  //
  if ( _dlist.size() > 0 ) {
    DetElemList::iterator diter = dlist().begin();
    DetElem* elem = 0;
    while( diter != dlist().end() ) {
      elem = *diter++;
      DetIntersection testinter = current;
      int flag = elem->intersect(traj,testinter);
      if ( flag && (testinter > current ) && (testinter < mininter) ) {
	mininter = testinter;
	//  reset the upper limit
	current.pathrange[1] = mininter.pathlen-_epsilon; 
	intersected = true;
      }
    }
  }
  //
  //  Apply the same recursively to all the subsets
  //
  if ( _slist.size() > 0 ) {
    DetSetList::iterator siter = slist().begin();
    DetSet* set = 0;
    while( siter != slist().end() ) {
      set = *siter++;
      DetIntersection setinter;
      if ( set->nextIntersection(traj,setinter,lastinter) && 
	   setinter < mininter ) {
	mininter = setinter;
	intersected = true;
      }
    }
  }
  nextinter = mininter;
  return intersected;
}
//
//  Another stupid navigation function, but more efficient than the above
//  for finding the full set of intersected elements
//

    
 void 
DchGasSet::intersection( std::vector<DetIntersection>& divec,
			 const Trajectory* traj,
			 double* myrange, bool clear ) const 
{
 //
 //  Make sure the set is ready
 //
 if ( !isReady() ) makeReady();
 //
 //  clear the vector
 //
 if ( clear ) divec.clear();
 //
 //  Find the range
 //
 double range[2];
 if ( myrange != 0 ) {
   range[0] = myrange[0];
   range[1] = myrange[1];
 } else {
   range[0] = traj->lowRange();
   range[1] = traj->hiRange();
 }
 //
 // cout<<"RANGE\t"<<range[0]<<"\t"<<range[1]<<endl;
 //
 //  Loop over the contained elements
 //
 if ( _dlist.size() > 0 ) {
   DetElemList::iterator diter = dlist().begin();
   DetElem* elem = 0;
   while( diter != dlist().end() ) {
     elem = *diter++;
     DetIntersection dinter(0,traj,range[0],range[1]);
     //
     //  Loop till the element is no longer intersected, to handle multiple 
     //  intersections
     //
     int ninter = 0;
     while ( elem->intersect(traj,dinter) && ninter < 2 ) {
	ninter++;
	if ( dinter.pathrange[0] < dinter.pathrange[1] ) {
	  divec.push_back(dinter);
	  dinter.pathrange[0] = dinter.pathrange[1] + _epsilon;
	  dinter.pathrange[1] = max(dinter.pathrange[0] + _epsilon,range[1]);
	} else {
	  break;
	}
	//
	//  Check if we're outside the range
	//
	if ( dinter.pathrange[0] > range[1] ) break;
     }
   }
 }
 std::sort(divec.begin(),divec.end());
 //
 //  Apply the same recursively to all the subsets without clearing
 //
 if ( _slist.size() > 0 ) {
   DetSetList::iterator siter = slist().begin();
   DetSet* set = 0;
   while( siter != slist().end() ) {
     set = *siter++;
     set->intersection(divec,traj,myrange,false);
   }
 }
}
    
