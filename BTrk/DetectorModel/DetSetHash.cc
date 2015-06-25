// ---------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetSetHash.cc,v 1.10 2004/12/14 07:10:18 bartoldu Exp $
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
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/ErrLog.hh"
#include "BTrk/DetectorModel/DetSetHash.hh"
#include "BTrk/DetectorModel/DetSet.hh"
#include <assert.h>
#include <iostream>

static const int _MAXINDEX = 1000000000; // arbitary large value

DetSetHash::DetSetHash(const DetSet& theSet,Ehash elemHash,elemSelFun esel,void* eseldat) :
  _hashtable(0),_maxindex(-_MAXINDEX),_minindex(_MAXINDEX),_nindex(0),_ehash(elemHash)
{
//
//  Make sure the set is locked
//
  assert(theSet.isLocked());
//
//  Make a list of all elements which pass the selection function
//
  DetElemList elemlist;
  bool clear= true;
  theSet.select(elemlist,esel,eseldat,0,0,clear);
//
//  Find the limits of the hash table
//
  DetElemList::const_iterator iter = elemlist.begin();
  DetElem* elem = 0;
  while ( iter != elemlist.end() ) {
    elem = *iter++;
    int elemnum = (*_ehash)(elem->elementNumber());
    _maxindex = std::max(_maxindex,elemnum);
    _minindex = std::min(_minindex,elemnum);
  }
//
//  Build the table
//
  _nindex = _maxindex-_minindex+1;
  _hashtable = new DetElem*[_nindex];
//
//  Initialize the array to zero
//
  for(int index=0;index<_nindex;index++)
    _hashtable[index] = 0;
//
//  Fill the table
//
  iter = elemlist.begin();
  while ( iter != elemlist.end() ) {
    elem = *iter++;
    int elemnum = (*_ehash)(elem->elementNumber());
//
//  Require at most one entry per cell
//
    assert(0 == _hashtable[elemnum-_minindex]);
    _hashtable[elemnum-_minindex] = elem;
  }
}
//
//  Access
//
DetElem*
DetSetHash::findElement(int elemidnum){
  int index = (*_ehash)(elemidnum);
  if(index >= _minindex && index <= _maxindex &&
     _hashtable[index] != 0 &&
     _hashtable[index]->elementNumber() == elemidnum)
    return _hashtable[index];
  else {
    ErrMsg(error) << "DetSetHash::findElement invalid element ID number" << endmsg;
    return 0;
  }
}
