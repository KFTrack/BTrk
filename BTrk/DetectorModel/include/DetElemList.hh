//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetElemList.hh,v 1.1 2004/12/14 07:10:17 bartoldu Exp $
//
// Description: 
//      Container for detector elements.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Rainer Bartoldus
//
// Copyright Information:
//      Copyright (C) 2004      Stanford Linear Accelerator Center
//
//--------------------------------------------------------------------------
#ifndef DETELEMLIST_HH
#define DETELEMLIST_HH

//---------------
// C++ Headers --
//---------------
//#include <list>
#include <vector>

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class DetElem;

//typedef std::list<DetElem*> DetElemList;
typedef std::vector<DetElem*> DetElemList;

#endif // DETELEMLIST_HH
