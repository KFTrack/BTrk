//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchMCGenTdcTime.hh 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchMCGenTdcTime
//      Base class for xxxxxxx
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Mirna van Hoek 
//	
//
// Copyright Information:
//	Copyright (C) 2001	University of Colorado
//
//------------------------------------------------------------------------

#ifndef DCHMCGENTDCTIME_HH
#define DCHMCGENTDCTIME_HH

class DchTimeToDist;

class DchMCGenTdcTime {
public:

 DchMCGenTdcTime();
 virtual  ~DchMCGenTdcTime();

 // raison d'etre...
 virtual double tdcTime(double gtime, double dist, int ambig, double z,
                 double charge, const DchTimeToDist& calib) const = 0;
 virtual const char* name() const = 0;

protected:
 bool genTime(double& dtime, double& dtimeRes,
                double gtime, double dist, 
                int ambig, double z, 
                const DchTimeToDist& calib) const;

private:

  //Note: if your class needs a copy constructor or an assignment operator, 
  //  make one of the following public and implement it.
  DchMCGenTdcTime( const DchMCGenTdcTime& );       // Copy Constructor
  DchMCGenTdcTime&       operator= ( const DchMCGenTdcTime& );  // Assignment op
};
#endif
