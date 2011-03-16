//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchWireStatus.cc 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchWireStatus
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Gerhard Raven           11/01/99
//
// Copyright Information:
//      Copyright (C) 1999      University of California, San Diego
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
//-----------------------
// This Class's Header --
//-----------------------
#include "DchCalib/DchWireStatus.hh"
#include <assert.h>
#include "ErrLogger/ErrLog.hh"


static const char rscid[] = "$Id: DchWireStatus.cc 88 2010-01-14 12:32:57Z stroili $";

DchWireStatus::DchWireStatus( unsigned layer, unsigned wire,
                              bool dead, bool noisy, bool hasBadCharge,
                              double eff, double effErr,
                              double noise, double noiseErr)
       : _layer(layer),_wire(wire),
         _eff(eff),_effErr(effErr),
         _noise(noise),_noiseErr(noiseErr),
         _noisy(noisy),_dead(dead),_badCharge(hasBadCharge)
{
        assert(layer>0);
}


DchWireStatus::~DchWireStatus()
{
}

double
DchWireStatus::efficiency() const
{
        if (_eff<0) 
             ErrMsg(fatal) << "To reduce CDB access, see http://babar-hn.slac.stanford.edu:5090/HyperNews/get/BdbSOS/1366/1/2/3/2.html, access to this data is suspended for now; a new interface will be written eventually" << endmsg;
        return _eff;
}

double
DchWireStatus::efficiencyErr() const
{ 
        if (_effErr<0) 
             ErrMsg(fatal) << "To reduce CDB access, see http://babar-hn.slac.stanford.edu:5090/HyperNews/get/BdbSOS/1366/1/2/3/2.html, access to this data is suspended for now; a new interface will be written eventually" << endmsg;
        return _effErr;
}

double
DchWireStatus::noise() const
{
        if (_noise<0) 
             ErrMsg(fatal) << "To reduce CDB access, see http://babar-hn.slac.stanford.edu:5090/HyperNews/get/BdbSOS/1366/1/2/3/2.html, access to this data is suspended for now; a new interface will be written eventually" << endmsg;
        return _noise;
}

double
DchWireStatus::noiseErr() const
{
        if (_noiseErr<0) 
             ErrMsg(fatal) << "To reduce CDB access, see http://babar-hn.slac.stanford.edu:5090/HyperNews/get/BdbSOS/1366/1/2/3/2.html, access to this data is suspended for now; a new interface will be written eventually" << endmsg;
        return _noiseErr;
}
