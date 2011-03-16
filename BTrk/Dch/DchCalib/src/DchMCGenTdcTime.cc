//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchMCGenTdcTime.cc 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchMCGenTdcTime
//      Base Class for
//      
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	M. van Hoek             
//	
//
// Copyright Information:
//	Copyright (C) 2001	University of Colorado
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include <math.h>

#include "DchCalib/DchMCGenTdcTime.hh"
#include "DchCalib/DchTimeToDist.hh"
#include "ErrLogger/ErrLog.hh"


DchMCGenTdcTime::DchMCGenTdcTime()
{ }

DchMCGenTdcTime::~DchMCGenTdcTime()
{}

bool 
DchMCGenTdcTime::genTime(double& dtime, double& dtimeRes,
                         double gtime, double dist,
                         int ambig, double z,
                         const DchTimeToDist &calib) const
{

  double x = calib.distToTime(dist,ambig,0,0,z,0);
  double eps = 0.001; //10 micron
  double alpha = (calib.distToTime(dist+eps,ambig,0,0,z,0) - x)/eps;
  if (alpha<0) return false; 
  dtime = x + gtime*1.e9; // gtime in s, dtime in ns!
  dtimeRes= alpha*calib.resolution(fabs(dist),ambig);
  return true;
}
