//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: ConsistencySet.hh 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//      A set of consistency objects. Typically used to define the genealogy 
//      of a CombinedConsistency object, i.e. contains a list of 
//      Consistencies the CombinedConsistency object is built from. 
//      The class provides methods to get an overlap, union etc. of two sets
//
//      The current implementation requires a string label to be attached 
//      to each Consistency object; it typically represents a PID subsystem 
//      "svt", "dch", etc.; case-insensitive.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Yury Kolomensky        05/02/98   
//          who shamelessly copied most of it from CombinedConsistency
//
//	Copyright (C) 1998      Caltech
//
//------------------------------------------------------------------------
#ifndef CONSISTENCYSET_HH
#define CONSISTENCYSET_HH

//-----------------
// BaBar Headers --
//-----------------
#include "BaBar/BaBar.hh"

//-------------
// C Headers --
//-------------
extern "C" {
#include <stddef.h>
}

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include <vector>
#include <string>

#include "ProbTools/Consistency.hh"

#include <iosfwd>

class ConsistencySet {
public:

  // default constructor; creates an empty set
  ConsistencySet();

  // copy constructor
  ConsistencySet(const ConsistencySet& rhs);

  virtual ~ConsistencySet();

  // assignment
  ConsistencySet& operator= (const ConsistencySet& rhs);

  // equality -- only checks the labels of parent objects, not their values
  virtual bool operator==(const ConsistencySet& rhs) const;
  virtual bool operator!=(const ConsistencySet& rhs) const;

  // modifiers
  
  // returns false if the operation cannot be performed (e.g., the label 
  // is already in the list
  virtual bool add( const char*, const Consistency& );

  // merge two sets. Tests for overlaps
  virtual bool combine(const ConsistencySet& other);

  // clear the lists
  virtual void reset();

  // accessors
  size_t nParents() const; 
  const Consistency* getConsistency( size_t index ) const;
  const Consistency* getLabelConsistency( const char* ) const;
  const char* getLabel( size_t index ) const;
  double      worstSignificance() const; 
  const char* worstSignificanceLabel() const;

  // overlaps -- returns a subset of self that overlaps with other
  // NOTE: the returned object is new'ed, client is responsible for deletion
  virtual ConsistencySet* overlap(const ConsistencySet& other) const;

  // print method
  virtual void print(std::ostream& os) const;

protected:

  //  helper functions
  bool   getLabelIndex_( const std::string&, size_t& index ) const;
  bool   worstSignificanceIndex_( size_t& index ) const;

private:
  std::vector<Consistency> _consistencyList;
  std::vector<std::string> _labelList;

};

#endif




