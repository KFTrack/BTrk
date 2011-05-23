//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchLayout.hh 91 2010-01-14 12:37:23Z stroili $
//
// Description:
//      DchLayout Class - pure virtual baseclass which defines
//      methods to obtain the general layout of the drift chamber
//      -- # of layers, and # of wires for each layer -- and 
//      no more than that! A lot of code doesn't need to know
//      more than that, and can thus depend on DchLayout instead
//      of on the full DchDetector (which should/will derive from DchLayout)
//
//      note: should add a method to obtain the # of the first wire,
//            and the # of the first layer
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Gerhard Raven
//
// Copyright Information:
//      Copyright (C) 2004              NIKHEF Amsterdam
//      
//
//------------------------------------------------------------------------

#ifndef DCHLAYOUT_HH
#define DCHLAYOUT_HH
#include <stddef.h>
class DchLayout {
public:
        virtual ~DchLayout() {;}
        virtual size_t nLayers() const =0;
        virtual size_t nWires(unsigned layer) const =0;
};
#endif
