//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchMCXTalk.hh 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchMCXTalk
//      Base class for specific X-talk models for DCH simulation
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	G. Raven        originator
//	
//
// Copyright Information:
//	Copyright (C) 2004	NIKHEF
//
//------------------------------------------------------------------------

#ifndef DCHMCXTALK_HH
#define DCHMCXTALK_HH


#include <string>
class DchDigi;
#include "CLHEP/Alist/AList.h"

class DchMCXTalk {
public:

 DchMCXTalk() {};
 virtual  ~DchMCXTalk() {};

 // raison d'etre...
 virtual void generateXTalk(HepAList<DchDigi>& digis) const = 0;
 virtual const std::string& name() const = 0;

protected:
 // utility function for subclasses to use
 // DchDigi *mergeXTalk(DchDigi* original, DchDigi* xtalk) const;

private:

  //Note: if your class needs a copy constructor or an assignment operator, 
  //  make one of the following public and implement it.
  DchMCXTalk( const DchMCXTalk& );       // Copy Constructor
  DchMCXTalk&       operator= ( const DchMCXTalk& );  // Assignment op
};
#endif // DCHMCXTALK_HH
