//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: CombinedConsistency.hh 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Gautier Hamel de Monchenault
//      Alexandre Telnov, December 2007: 
//         add logLikelihood 
//
// Copyright Information:
//	Copyright (C) 1998, 2007
//
//------------------------------------------------------------------------
#ifndef COMBINEDCONSISTENCY_HH
#define COMBINEDCONSISTENCY_HH

//---------------
// C++ Headers --
//---------------

//----------------------
// Base Class Headers --
//----------------------
#include "ProbTools/Consistency.hh"
class ConsistencySet;

class CombinedConsistency  : public Consistency 
{
  
public:

  // constructors from the genealogy  object

  // this one takes over ownership
  CombinedConsistency(ConsistencySet* genealogy);

  // and this one makes own copy
  CombinedConsistency(const ConsistencySet& genealogy);

  // copy constructor; 
  CombinedConsistency( const CombinedConsistency& );

  // assignement operator
  CombinedConsistency& operator=( const CombinedConsistency& );

  // destructor
  virtual ~CombinedConsistency();

  // genealogy (base class virtual functions)
  virtual const ConsistencySet* genealogy() const {return _genealogy;}

  // print method
  virtual void print(std::ostream& os) const;

protected:
  // default ctor -- to be used only by subclasses that
  // override combining of consistencies
  CombinedConsistency();

private:
  // compute the combined likelihood
  void combine();
  
private:
  
  ConsistencySet* _genealogy;
  
};

#endif




