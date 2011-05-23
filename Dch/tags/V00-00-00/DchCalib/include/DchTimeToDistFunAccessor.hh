#ifndef DCHTIMETODISTFUNACCESSOR_HH
#define DCHTIMETODISTFUNACCESSOR_HH

//--------------------------------------------------------------------------
//
// Environment:
//      This software was developed for the BaBar collaboration.  If you
//      use all or part of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 2000      University of California, San Diego
//
//------------------------------------------------------------------------


/**
 *  DchTimeToDistFunAccessor. 
 *     this class regulates access to some protected functions
 *     in DchTimeToDist: only classes inheriting from this one can 
 *     access the actual functions used for t->d calibrations
 *
 *
 *  This software was developed for the BaBar collaboration.  If you
 *  use all or part of it, please give an appropriate acknowledgement.
 *
 *  Copyright (C) 2000 University of California, San Diego
 *
 *  @version $Id: DchTimeToDistFunAccessor.hh 88 2010-01-14 12:32:57Z stroili $
 *
 *  @author (Gerhard Raven)           (based on an idea of Steve Schaffner)
 */

#include "DchCalib/DchTimeToDist.hh"

class DchTimeToDistFunAccessor
{
public:
  virtual ~DchTimeToDistFunAccessor() = 0;
protected:
  double t0(const DchTimeToDist&x) const { return x.t0();}
  const DchTimeToDistFun* timeToDistFun(const DchTimeToDist&x) const
         { return x.timeToDistFun();}
  const DchCalibFun* resolutionFun(const DchTimeToDist&x) const
         { return x.resolutionFun();}
};

#endif
