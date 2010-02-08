// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetAEAdder.hh,v 1.1 2002/01/15 21:34:23 brownd Exp $
//
//  Description:  DetAEAdder. Implementation of DetAECombiner that 'adds' two
//  alignments.  Note that the covariance matrix is taken from the second (delta)
//  set if alignment parameters
//
// Copyright Information:
//	Copyright (C) 2002	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 1/15/2002
//------------------------------------------------------------------------------
#ifndef DETAEADER_HH
#define DETAEADER_HH

#include "DetectorModel/DetAECombiner.hh"

class DetAEAdder : public DetAECombiner {
public:
// default constructor only
  DetAEAdder() {}
  ~DetAEAdder() {}
// only one function: combine two DetAlignElems
  virtual bool combine(const DetAlignElem& inone,
		       const DetAlignElem& intwo,
		       DetAlignElem& out) const;
private:
// preempt copy and equivalence
  DetAEAdder& operator = (const DetAEAdder&);
  DetAEAdder(const DetAEAdder&);
};


#endif
