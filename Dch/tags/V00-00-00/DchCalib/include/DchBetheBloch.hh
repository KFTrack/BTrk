#ifndef DCHBETHEBLOCH_HH
#define DCHBETHEBLOCH_HH
//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchBetheBloch.hh 132 2010-09-21 11:15:12Z stroili $
//
// Description:
//	Class DchBetheBloch
//	Dch Transient class containing DchPID Bethe-Bloch parameters
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	F. Galeazzi		Originator
//
// Copyright Information:
//	Copyright (C) 1999	INFN & University of Padova
//
// Revision History:
//	20070328  Move ::getDchMomentum(TrkFit*) here from DchPidInfo
//	20070410  Restore original ::ionization() interface, adding TRT*
//		  as third arg, defaulting to null.
//------------------------------------------------------------------------

#include "PDT/PdtPid.hh"
#include <vector>
#include <string>
#include <iosfwd>

class TrkRecoTrk;
class TrkFit;

class DchBetheBloch 
{
public:
  // Constructors
  DchBetheBloch( const std::vector<double>& parVector,
                 const char* functionName );
  // Destructor
  virtual ~DchBetheBloch( );

  const char* channelName( ) const { return _functionName.c_str(); }
  unsigned numberOfParameters( ) const { return _par.size(); }
  double fieldValue(unsigned i) const { return _par[i]; }

  virtual double ionization(double momentum, PdtPid::PidType hypo,
			    const TrkRecoTrk* track=0) const = 0;

  void print( std::ostream& o) const;

private:
  // Data members
  std::vector<double> _par;
  std::string _functionName;
  DchBetheBloch();
};

#endif // DCHBETHEBLOCH_HH
