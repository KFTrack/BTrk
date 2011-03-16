// ----------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchSelections.cc 91 2010-01-14 12:37:23Z stroili $
//
//  Description:
//  This contains functions for making various selections of Dch elements
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown and Gerry Lynch
// ----------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include <assert.h>

#include "DchGeom/DchSelections.hh" 
#include "DchGeom/DchSearchId.hh"
#include "DetectorModel/DetSet.hh"
#include <string>

//
//  ID selection functions
//

// select Inner and Outer Cylinder
bool selectDchInner(const DetElem* elem,void* data){
  return (elem->elementName().find("Dch Inner") != std::string::npos );
}

bool selectDchOuter(const DetElem* elem,void* data){
  return (elem->elementName().find("Dch Outer") != std::string::npos );
}

// select rear end plate
bool selectDchRPlate(const DetElem* elem,void* data){
  return (elem->elementName().find("Dch rear") != std::string::npos );
}

// select forward End Plate
bool selectDchFPlate(const DetElem* elem,void* data){
  return (elem->elementName().find("Dch forward") != std::string::npos );
}

// select gas volume element
bool selectDchGasVol(const DetElem* elem,void* data){
  return (elem->elementName().find("Dch volume gas") != std::string::npos );
}

bool selectDchSetID(const DetSet* set,void* data){
  DchSearchId* sid = (DchSearchId*)data;
  if(sid)
    return (set->setName().find(sid->tag->c_str()) != std::string::npos ) &&
    set->setNumber() >= sid->idLim[0] &&
    set->setNumber() <= sid->idLim[1];
  else
    return true;
}

// find layer elements
bool findDchLayer(const DetElem* elem,void* data){
  return (elem->elementName().find("Dch Layer") != std::string::npos ); // make sure it's a layer!
}

// find SuperLayer elements
bool findDchSuperLayer(const DetSet* set, void* data) {
  return (bool) (set->setName().find("Dch SuperLayer") != std::string::npos ); // make sure it's a SuperLayer!
}

// find any Dch element
bool selectDchElement(const DetElem* elem, void* data){
  return (elem->elementName().find("Dch") != std::string::npos );
}

bool selectDchSet(const DetSet* set,void* data){
  DchSearchId* sid = (DchSearchId*)data;
  return   ( set->setName().find(sid->tag->c_str()) != std::string::npos)   &&
           set->setNumber() >= sid->idLim[0]          &&
           set->setNumber() <= sid->idLim[1];
}
