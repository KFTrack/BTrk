// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetAECombiner.hh,v 1.1 2002/01/15 21:34:23 brownd Exp $
//
//  Description:  DetAECombiner.  An abstract factory for combining 2 DetAlignElem
//  objects to create a third.
//
// Copyright Information:
//	Copyright (C) 2002	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 1/15/2002
//------------------------------------------------------------------------------
#ifndef DETAECOMBINER_HH
#define DETAECOMBINER_HH

class DetAlignElem;

class DetAECombiner {
public:
// default constructor only
  DetAECombiner() {}
  virtual ~DetAECombiner() {}
// only one function: combine two DetAlignElems
  virtual bool combine(const DetAlignElem& inone,
		       const DetAlignElem& intwo,
		       DetAlignElem& out) const = 0;
private:
// preempt copy and equivalence
  DetAECombiner& operator = (const DetAECombiner&);
  DetAECombiner(const DetAECombiner&);
};


#endif
