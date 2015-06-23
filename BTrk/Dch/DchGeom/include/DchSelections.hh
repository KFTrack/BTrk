// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchSelections.hh 91 2010-01-14 12:37:23Z stroili $
//
//  Description:
//  This contains functions for making various selections of Dch elements
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown and Gerry Lynch
// ------------------------------------------------------------------------------
#ifndef DCHSELECTIONS_HH
#define DCHSELECTIONS_HH

#include <iostream>
#include <fstream>

class DetElem;
class DetSet;
//
//  Selection functions for finding Svt elements
//
bool selectDchInner(const DetElem* elem,void* data);
bool selectDchOuter(const DetElem* elem,void* data);
bool selectDchRPlate(const DetElem* elem,void* data);
bool selectDchFPlate(const DetElem* elem,void* data);
bool selectDchGasVol(const DetElem* elem,void* data);
bool selectDchSet(const DetSet* set,void* data);
bool selectDchSetID(const DetSet*,void*);
bool findDchLayer(const DetElem*,void*);
bool findDchSuperLayer(const DetSet* set, void*);
bool selectDchElement(const DetElem*,void*);

#endif
