// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetExternalSet.hh,v 1.6 2002/12/17 04:03:57 dbrown Exp $
//
//  Description:
//  Special form of a DetSet where the elements are outside the
//  tracking volume.  For this set the intersection methods are trivial.
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
//  Authors: Gautier Hamel de Monchenault, 1/7/97
//------------------------------------------------------------------------------
#ifndef EXTERNALSET_HH
#define EXTERNALSET_HH
//
//  Includes
//
#include "BTrk/DetectorModel/DetSet.hh"
#include <vector>

//
//  Define the class
//
class DetExternalSet : public DetSet {
public:
//  Unique constructor
  DetExternalSet(const char* name,int IDNum);
//  Destructor
  virtual ~DetExternalSet() {;}
//
//  Just overwrite the intersection methods
//
  bool firstIntersection(const Trajectory*,DetIntersection& next,
			 double* myrange = 0) const;
  bool nextIntersection(const Trajectory*,DetIntersection& next,
                        DetIntersection* prev = 0) const;
  void intersection(std::vector<DetIntersection>&,
		    const Trajectory*,double* myrange=0,
		    bool clear=true) const;

};
#endif
