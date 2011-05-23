//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchWireStatusList.cc 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchWireStatusList
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	G. Raven                
//	
//
// Copyright Information:
//	Copyright (C) 1999	University of California, San Diego
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include "DchCalib/DchWireStatusList.hh"
#include "DchCalib/DchWireStatus.hh"
#include "ErrLogger/ErrLog.hh"

DchWireStatusList::DchWireStatusList(std::auto_ptr<DchWireArray<DchWireStatusRCPtr > >& status)
        :_status(status.release())
{
}

DchWireStatusList::~DchWireStatusList()
{
}
