// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetSet.hh,v 1.25 2004/12/14 07:10:18 bartoldu Exp $
//
//  Description:
//  Base Class to define a (heirarchical) list of DetElems.  A DetSet
//  object can contain (either,or, or both) a list of DetElems and
//  a list of other DetSets.  These sets provide the navigation tools for
//  the needs of tracking.
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 10/10/96
//------------------------------------------------------------------------------
#ifndef DETECTORSET_HH
#define DETECTORSET_HH
//----------------
// BaBar header --
//----------------
#if defined( HP1022 ) && !defined( BABAR_HH )
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/ErrLog.hh"
#endif // HP1022 && !BABAR_HH
//
//  global includes
//
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
//
//  BaBar includes
//
#include "BTrk/DetectorModel/DetIntersection.hh"
#include "BTrk/DetectorModel/DetElemList.hh"
#include "BTrk/DetectorModel/DetSetList.hh"
class DetSet;
class DetElem;
class Trajectory;
class DetAlignElem;
//
//  List function prototypes
//
typedef bool (*elemSelFun)(const DetElem*,void*);
typedef bool (*setSelFun)(const DetSet*,void*);
//
void dsElemPrint(const DetElem*,void*); // simple printout
void dsSetPrintAll(const DetSet*,void*);
void dsGnuplotElement(const DetElem*,void*); // gnuplot the elements
bool dsMatchAll(const DetElem*,void*); // find all elements
//
//  Define the main class
//
class DetSet{
public:
//
//  Constructors
//
  DetSet(const char*,int);
  DetSet();
//  Destructor
  virtual ~DetSet();
//
//  Operators
//
  DetSet& operator += (DetElem& elem);
  DetSet& operator += (DetSet& set);
//  equality
  bool operator == (const DetSet& other){
    return _dsnum == other._dsnum && _dsname == other._dsname; }
//
//  Access
//
  virtual void print(std::ostream& os) const;
  virtual void printAll(std::ostream& os ) const;
  int setNumber() const { return _dsnum;}
  const std::string& setName() const { return _dsname;}
//
//  Alignment functions. Global alignment is 
//  passed down recursively to all contained elements.  
//  Local alignment is applied to all contained detectorelements
//  which match the input alignElements.  The return value
//  tells how many input DetAlignElems found a matching
//  DetElem. These functions
//  are virtual to allow derrived classes to update their local
//  data members by overwriting this function, but a pedantic implementation
//  is provided.
//
  virtual void applyGlobal(const DetAlignElem&);// apply global alignment
  virtual void removeGlobal(const DetAlignElem&);// unapply global alignment
// please use the following new alignment functions:
  virtual int applyLocal(const std::vector<DetAlignElem>&);
  virtual int removeLocal(const std::vector<DetAlignElem>&);
// the following are deprecated
  virtual int applyLocal(const DetAlignElem*,int nelem);
  virtual int removeLocal(const DetAlignElem*,int nelem);
//
//  Recursively select elements according to input functions, returning
//  a list of selected element pointers.  The user must supply the
//  element selection function, but the subset selection
//  function is optional (it can be provided to improve efficiency for
//  nested elements).  The last argument controls whether the input list
//  is zeroed, and should normally be left to its default value
//  (it's used in recursive sub-calls).
//  
  void select
  ( DetElemList&, elemSelFun, void* elemdata,
    setSelFun setfun = 0,void* setdata= 0,bool clear = true) const;
//
//  Same for entire sets
//
  void setSelect
  ( DetSetList&, setSelFun, void* setdata,
    bool clear = true) const;
//  Access to directly contained elements/sets.
  const DetElemList& elementList() const { return _dlist; }
  const DetSetList&      setList() const { return _slist; }
//
//  Find the next element intersected by a given Trajectory from a given starting intersection
//  Don't use this to loop over elements, as it is inefficient.
//  The default implementation of this calls down to firstIntersection,
//  so subclasses don't need to override this.
//
  virtual bool nextIntersection(const Trajectory* traj,DetIntersection& next,
				DetIntersection* previous = 0) const;
//
//  Similiar function, which finds the first intersection within the
//  trajectory range (or the explicitly specified range).  This is the
//  function smarter subclasses of DetSet should override.
//
  virtual bool firstIntersection(const Trajectory* traj,DetIntersection& next,
				 double* myrange=0) const;
//
//  Make a complete list of interesected DetElems for a given trajectory.  The (optional)
//  supplied range overrides the trajectory range
//
  virtual void intersection(std::vector<DetIntersection>&,
			    const Trajectory*,double* myrange=0,
			    bool clear=true) const;

//
  virtual void setLock(); // Recursively prevent new elements or sets from being added.
  bool isLocked() const { return _lock; } // access to lock status
  virtual bool isReady() const; // check readyness of a set
  virtual void makeReady() const; // Make a set ready

protected:
  int _dsnum; // define the set number
  std::string _dsname; // define the set name
  DetElemList _dlist; // pointer list to elements
  DetSetList  _slist; // pointer list to (recursive) lists
  bool _lock; // flag to lock the set
  bool _ready; // flag to define the set as ready for geometric operations

public:
// Returns a pointer list of the elements of this set (if clear=true)
// or appends to the list the elements of this set (if clear=false)
// The list includes all the elements of subsets, recursively
  void listAllElements(DetElemList& list, 
		       bool clear=true) const {
    select(list,dsMatchAll,0,0,0,clear); }

protected:
//
//  Functions to cast of 'const'; this should be provided by the language
//
  DetElemList& dlist() const {
    return (DetElemList&)_dlist;}
  DetSetList& slist() const {
    return (DetSetList&)_slist; }
  bool& lock() const { return (bool&) _lock;}
  bool& ready() const { return (bool&) _ready; }
};

#endif
