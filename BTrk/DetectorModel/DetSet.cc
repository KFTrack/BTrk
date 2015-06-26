// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetSet.cc,v 1.38 2004/12/14 07:10:18 bartoldu Exp $
//
//  Description:  Class DetSet
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 10/10/96
//------------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"
#include <vector>
#include "BTrk/BaBar/ErrLog.hh"
#include <assert.h>
#include "BTrk/DetectorModel/DetIntersection.hh"
#include "BTrk/DetectorModel/DetSet.hh"
#include "BTrk/DetectorModel/DetElem.hh"
#include "BTrk/DetectorModel/DetAlignElem.hh"
#include "BTrk/BbrGeom/Trajectory.hh"
using std::endl;
using std::ostream;
static const double _epsilon = 1.0e-4;
static const int _MAXINTER = 3000; // maximum # of successive intersections of a single element
//
//  Define the list functions for printing
//
void dsElemPrint(const DetElem* elem,void* os){
  ostream* stream = (ostream*)os;
  elem->print(*stream);
}

void dsSetPrintAll(const DetSet* set,void* os){
  ostream* stream = (ostream*)os;
  set->printAll(*stream);
}

//  Id selection function (not a member)
bool dsSelIDFun(const DetElem* elem,void* data){
  return elem->elementNumber() == *(int*)data;
}
//
//  Trivial match function
//
bool dsMatchAll(const DetElem*,void*){
  return true; }

//
//  Define the class functions
//
DetSet::DetSet(const char* name,int iset)
  : _dsnum(iset), _dsname(name), _lock(false), _ready(false)
{
}

DetSet::DetSet()
  : _dsnum(-1), _dsname("Unknown"), _lock(false), _ready(false)
{
}

//
DetSet::~DetSet(){
}
//
//  Operators
//
DetSet&
DetSet::operator+=(DetElem& elem) {
  assert(!_lock);
  _dlist.push_back(&elem);
  return *this;
}
//
DetSet&
DetSet::operator+=(DetSet& set) {
  assert(!_lock);
  _slist.push_back(&set);
  return *this;
}
//
//  trivial implementations for testing readiness and setting the lock.
//
bool
DetSet::isReady() const {
  return _ready;
}
void
DetSet::makeReady() const {
  ready() = true;
//
//  Recursively make subsets ready
//
  if(_slist.size() > 0){
    DetSetList::const_iterator siter = slist().begin();
    DetSet* set = 0;
    while ( siter != slist().end() ) {
      set = *siter++;
      set->makeReady();
    }
  }
}

void
DetSet::setLock(){
  lock() = true;
//
//  Recursively set locks on the subsets
//
  if(_slist.size() > 0){
    DetSetList::const_iterator siter = slist().begin();
    DetSet* set = 0;
    while ( siter != slist().end() ) {
      set = *siter++;
      set->setLock();
    }
  }
}
//
//  Print; this uses two list functions to access the elements.
//
void
DetSet::print(ostream& os) const{
  os << "Detector Set " << _dsname << " " << _dsnum
     << " has " << _dlist.size() << " direct elements and "
     << _slist.size() << " direct subsets";
  if(isReady())
    os << " and is Ready for intersections" << endl;
  else
    os << " and is NOTReady for intersections" << endl;
}
//
//  Print; this uses two list functions to access the elements.
//
void
DetSet::printAll(ostream& os) const{
//
//  First, check for direct elements
//
  if(_dlist.size() > 0){
    os << "Detector Set " << _dsname << " " << _dsnum << 
      " has the following direct elements: "<< endl;
    DetElemList::const_iterator diter = dlist().begin();
    DetElem* anElem = 0;
    while ( diter != dlist().end() ) {
      anElem = *diter++;
      dsElemPrint( anElem, (void*)&os );
      //    dlist().apply(&dsElemPrint,(void*)&os);
    }
  } else
    os << "Detector Set " << _dsname << " " << _dsnum <<
      " has no direct elements ";
  if(_slist.size() > 0){
    os << "and contains the following subsets: "<< endl;
    DetSetList::const_iterator siter = slist().begin();
    DetSet* aSet = 0;
    while ( siter != slist().end() ) {
      aSet = *siter++;
      dsSetPrintAll( aSet, (void*)&os );
      //    slist().apply(&dsSetPrintAll,(void*)&os);
    }
  } else
    os << "and contains no subsets "<< endl;
}
//
//  Recurisively select elements according to user supplied functions.
//  The output is returned as a list of elements

//STL version

void
DetSet::select( DetElemList& theList,
		elemSelFun efun, void* elemdata,
		setSelFun sfun, void* setdata,
		bool clear ) const
{
  //
  //  clear the list at the top level call
  //
  if ( clear )
    theList.clear();
  //
  //  Select from the contained elements of this set if it is selected or there
  //  is no set selection function
  //
  if ( !sfun || (*sfun)( this, setdata ) )
    {
      for ( DetElemList::const_iterator diter = _dlist.begin();
	    diter != _dlist.end();
	    ++diter )
	{
	  DetElem* match = *diter;
	  if ( efun( match, elemdata ) )
	    theList.push_back( match );
	}
    }
  //
  //  Recursively select on the contained sets, WITHOUT clearing the list!
  //
  for ( DetSetList::const_iterator siter = _slist.begin();
	siter != _slist.end();
	++siter )
    {
      const DetSet* match = *siter;
      match->select( theList, efun, elemdata, sfun, setdata, false );
    }
}


//
//  Recursively select sets according to user supplied function.
//  The output is returned as a list of sets.
//
void
DetSet::setSelect( DetSetList& theList,
		   setSelFun sfun, void* setdata,
		   bool clear ) const
{
  //
  //  clear the list at the top level call
  //
  if ( clear )
    theList.clear();
  //
  //  Check this set
  //
  if ( (*sfun)( this, setdata ) )
    theList.push_back( (DetSet*)this );
  //
  //  Recursively select on the contained sets, WITHOUT clearing the list!
  //
  for ( DetSetList::const_iterator siter = _slist.begin();
	siter != _slist.end();
	++siter )
    {
      const DetSet* match = *siter;
      match->setSelect( theList, sfun, setdata, false );
    }
}


//
//  Navigation; find the 'next' element, given a trajectory and a current
//  path length.  This is the stupidest imaginable implementation, trying
//  to force an intersection with every element, and testing on the return
//  value.  Smarter sets should overwrite this function!!!
//
bool
DetSet::nextIntersection(const Trajectory* traj,
			      DetIntersection& nextinter,
			      DetIntersection* lastinter) const {
//
// Call down to firstIntersection with the appropriate arguments
//
  if(lastinter){
    double trajrange[2];
    trajrange[0] = lastinter->pathlen + _epsilon;
    trajrange[1] = traj->hiRange();
    return firstIntersection(traj,nextinter,trajrange);
  } else
    return firstIntersection(traj,nextinter);
}


bool
DetSet::firstIntersection(const Trajectory* traj,
			  DetIntersection& firstinter,
			  double* myrange)  const {
//
//  Make sure the set is ready
//
  if(!isReady())makeReady();
//
//  use the explicit range if provided,
//  else start from the begining of the trajectory range
//
  bool intersected = false;
  double trajrange[2];
  if(myrange != 0){
    trajrange[0] = myrange[0];
    trajrange[1] = myrange[1];
  } else {
    trajrange[0] = traj->lowRange();
    trajrange[1] = traj->hiRange();
  }
  DetIntersection mininter(0,traj,trajrange[1]);
//
//  Loop over the contained elements, and intersect the trajectory with them
//
  if(_dlist.size() > 0){
    DetIntersection current(0,traj,trajrange[0],trajrange[1]);
    DetElemList::const_iterator eiter = dlist().begin();
    DetElem* elem = 0;
    while ( eiter != dlist().end() ) {
      elem = *eiter++;
      DetIntersection testinter(current);
      if(elem->intersect(traj,testinter)) {
	mininter = testinter;
	current.pathrange[1] = mininter.pathlen-_epsilon; // reset the upper limit
	intersected = true;
// abort if we're at the limit; can't do better than this!
	if(current.pathrange[0] >= current.pathrange[1]){
	  firstinter = mininter;
	  return intersected;
	}
      }
    }
// reset the range
    if(intersected)
      trajrange[1] = current.pathrange[1];
  }
//
//  Apply the same recursively to all the subsets
//
  if(_slist.size() > 0){
    DetSetList::const_iterator siter = slist().begin();
    DetSet* set = 0;
    while ( siter != slist().end() ) {
      set = *siter++;
      DetIntersection setinter;
      if(set->firstIntersection(traj,setinter,trajrange)){
	mininter = setinter;
	intersected = true;
	trajrange[1] = mininter.pathlen-_epsilon;
      }
    }
  }
  firstinter = mininter;
  return intersected;
}
//
//  Another stupid navigation function, but more efficient than the above
//  for finding the full set of intersected elements
//

void
DetSet::intersection(std::vector<DetIntersection>& divec,
		     const Trajectory* traj,
		     double* myrange,bool clear) const {
//
//  Make sure the set is ready
//
  if(!isReady())makeReady();
//
//  clear the vector
//
  if(clear)divec.clear();
//
//  Find the range
//
  double range[2];
  if(myrange != 0){
    range[0] = myrange[0];
    range[1] = myrange[1];
  } else {
    range[0] = traj->lowRange();
    range[1] = traj->hiRange();
  }
//
//  Loop over the contained elements
//
  if(_dlist.size() > 0){
    DetElemList::const_iterator diter = dlist().begin();
    DetElem* elem = 0;
    while ( diter != dlist().end() ) {
 
      elem = *diter++;
      if (ErrLogging(debugging)) {
	ErrMsg(debugging) << "begin "<<elem->elementName()<<endmsg;
      }
      DetIntersection dinter(0,traj,range[0],range[1]);
//
//  Loop till the element is no longer intersected, to handle multiple intersections
//
      int ninter = 0;
      while(elem->intersect(traj,dinter) && ninter < _MAXINTER){
	ninter++;
// require a valid intersection
	if(dinter.pathrange[0] < dinter.pathrange[1]){
	  divec.push_back(dinter);
// abort if we're at the end of the range
	  if(dinter.pathrange[1] >= range[1] )
	    break;
// setup for the next intersectoin
  	  dinter.pathrange[0] = dinter.pathrange[1] + _epsilon;
          dinter.pathrange[1] = std::max(dinter.pathrange[0] + _epsilon,range[1]);
	} else {
	  ErrMsg(error) << "DetSet: exit point comes before entrance in element interection." << endmsg;
	  break;
	}
      }
    }
  }
  std::sort(divec.begin(), divec.end() );

//
//  Apply the same recursively to all the subsets without clearing
//
  if(_slist.size() > 0){
    DetSetList::const_iterator siter = slist().begin();
    DetSet* set = 0;
    while ( siter != slist().end() ) {
      set = *siter++;
      set->intersection(divec,traj,myrange,false);
    }
  }
}


//
//  Generic implementations of alignment functions
//
void
DetSet::applyGlobal(const DetAlignElem& glob){
//
//  Loop over the contained elements and apply the alignment
//
  if(_dlist.size() > 0){
    DetElemList::const_iterator diter = dlist().begin();
    DetElem* elem = 0;
    while ( diter != dlist().end() ) {
      elem = *diter++;
      elem->applyGlobal(glob);
    }
  }
//
//  Apply the same recursively to all the subsets
//
  if(_slist.size() > 0){
    DetSetList::const_iterator siter = slist().begin();
    DetSet* set = 0;
    while ( siter != slist().end() ) {
      set = *siter++;
      set->applyGlobal(glob);
    }
  }
//
//  After alignment, a set is no longer ready
//
  _ready = false;
}
//
void
DetSet::removeGlobal(const DetAlignElem& glob){
//
//  Loop over the contained elements and remove the alignment
//
  if(_dlist.size() > 0){
    DetElemList::const_iterator diter = dlist().begin();
    DetElem* elem = 0;
    while ( diter != dlist().end() ) {
      elem = *diter++;
      elem->removeGlobal(glob);
    }
  }
//
//  Apply the same recursively to all the subsets
//
  if(_slist.size() > 0){
    DetSetList::const_iterator siter = slist().begin();
    DetSet* set = 0;
    while ( siter != slist().end() ) {
      set = *siter++;
      set->removeGlobal(glob);
    }
  }
  _ready = false;
}
int
DetSet::applyLocal(const DetAlignElem* alist, unsigned nelem){
  ErrMsg(warning) << "This is an obsolete function: please convert to call " << endl
		  << "DetSet::applyLocal(std::vector<DetAlignElem>) instead" << endmsg;
  std::vector<DetAlignElem> elems;
  for(unsigned ielem=0;ielem<nelem;ielem++)
    elems.push_back(alist[ielem]);
  return applyLocal(elems);
}
//
//  Same thing for removeal
//
int
DetSet::removeLocal(const DetAlignElem* alist, unsigned nelem){
  ErrMsg(warning) << "This is an obsolete function: please convert to call " << endl
		  << "DetSet::removeLocal(std::vector<DetAlignElem>) instead" << endmsg;
  std::vector<DetAlignElem> elems;
  for(unsigned ielem=0;ielem<nelem;ielem++)
    elems.push_back(alist[ielem]);
  return removeLocal(elems);
}


int
DetSet::applyLocal(const std::vector<DetAlignElem>& alist){
  int nfound = 0;
//
//  Loop over the list of alignments, find the matching
//  detectorElement
//
  DetElemList::const_iterator diter = _dlist.begin();
  while(diter != _dlist.end()){
// try to find a match with the align elements for this element
    std::vector<DetAlignElem>::const_iterator ae = alist.begin();
    while(ae != alist.end()){
      if((*diter)->match(*ae)){
	nfound++;
	(*diter)->applyLocal(*ae);
	break;
      }
      ae++;
    }
    diter++;
  }
//
//  Apply recursively to the sublists
//
  DetSetList::const_iterator siter = _slist.begin();
  while(siter != _slist.end()){
    nfound += (*siter)->applyLocal(alist);
    siter++;
  }
// mark the cache invalid
  _ready = false;
  return nfound;
}
//
//  Same thing for removeal
//

int
DetSet::removeLocal(const std::vector<DetAlignElem>& alist){
  int nfound = 0;
//
//  Loop over the list of alignments, find the matching
//  detectorElement
//
  DetElemList::const_iterator diter = _dlist.begin();
  while(diter != _dlist.end()){
// try to find a match with the align elements for this element
    std::vector<DetAlignElem>::const_iterator ae = alist.begin();
    while(ae != alist.end()){
      if((*diter)->match(*ae)){
	nfound++;
	(*diter)->removeLocal(*ae);
	break;
      }
      ae++;
    }
    diter++;
  }
//
//  Apply recursively to the sublists
//
  DetSetList::const_iterator siter = _slist.begin();
  while(siter != _slist.end()){
    nfound += (*siter)->removeLocal(alist);
    siter++;
  }
// mark the cache invalid
  _ready = false;
  return nfound;
}

//
//  Function to gnuplot an element, in list applicable form
//
void dsGnuplotElement(const DetElem* elem,void* plot){
  GnuPlot* gp = (GnuPlot*) plot;  // cast plot as a GnuPlot
  elem->gnuPlot( gp );            // forward to the element
}
