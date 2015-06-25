// ---------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetSetHash.hh,v 1.4 1999/12/18 17:14:26 brownd Exp $
//
//  Description:  A hash table for DetSet objects.  This allows efficient
//  access to elements according to a (user supplied) hash function based
//  on their ID number.
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 3/1/97
//----------------------------------------------------------------------------

#ifndef DETSETHASH_HH
#define DETSETHASH_HH

#include "BTrk/DetectorModel/DetSet.hh"

typedef int (*Ehash)(int); // hash function prototype

class DetSetHash {
public:
//
//  Construct from a (locked) DetSet object, given the element ID# hashing
//  function.  This function has to return a unique integer for each
//  set ID number, and the output (over the set) should be as compact as
//  possible to avoid wasting memory
//
  DetSetHash(const DetSet&,Ehash,elemSelFun,void* data = 0);
// destructor
  ~DetSetHash(){ if(_hashtable != 0) delete[] _hashtable; }
//
//  Only one access function: find the element given it's ID number.
//
  DetElem* findElement(int elemidnum);
//
private:
//  Array of DetElems indexed by hash function applied to their
//  ID number
  DetElem** _hashtable;
  int _nindex;
  int _maxindex;
  int _minindex; // define the index space
  Ehash _ehash; // save a pointer to the hash function
// prohibit
  DetSetHash& operator = (const DetSetHash&);
  DetSetHash(const DetSetHash&);
};
#endif
