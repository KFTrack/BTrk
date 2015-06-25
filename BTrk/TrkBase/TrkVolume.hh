//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkVolume.hh,v 1.9 2003/01/21 12:55:11 raven Exp $
//
// Description:
//	TrkVolume Class -  
//      Abstract interface for the description
//      of the tracking volume of any subdetector
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Gautier Hamel de Monchenault - CEN Saclay & Lawrence Berkeley Lab
//
// History (add to end):
//      Gautier   May 6, 1997  - creation
//      Dave Brown 5/11/97  Add a name member
//
// Copyright Information:
//	Copyright (C) 1997		Lawrence Berkeley Laboratory
//	Copyright (C) 1997	       CEA - Centre d'Etude de Saclay
//
//------------------------------------------------------------------------

#ifndef TRKVOLUME_HH
#define TRKVOLUME_HH

#include "BTrk/TrkBase/TrkDirection.hh"
#include <string>

class Trajectory;
class HepPoint;

class TrkVolume {

public:
// Default ctor
  TrkVolume();
  TrkVolume(const char*);
// Destructor
  virtual ~TrkVolume();
//
  virtual bool extendThrough( const Trajectory* theTraj,
			      double& theFlightDist,
			      trkDirection theDirection=trkOut,
			      double* theStartingFlightDist=0) const = 0;
// returns false if the point is outside the volume
  virtual bool isInside( const HepPoint& ) const = 0;
// Access to the name
  const std::string& name() const {return _tvname;}
private:
  std::string _tvname; // define the volume name
};

#endif
