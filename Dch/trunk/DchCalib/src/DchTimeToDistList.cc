//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchTimeToDistList.cc 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchTimeToDistList
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	G. Raven                
//	
//
// Copyright Information:
//	Copyright (C) 1998	University of California, San Diego
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include "DchCalib/DchTimeToDistList.hh"
#include "DchCalib/DchTimeToDist.hh"

static const char rcsid[] = "$Id: DchTimeToDistList.cc 88 2010-01-14 12:32:57Z stroili $";

DchTimeToDistList::DchTimeToDistList(std::auto_ptr<DchWireArray<DchTimeToDistRCPtr> >& t2d)
                :_t2d(t2d.release())
{
}

DchTimeToDistList::~DchTimeToDistList()
{
}
