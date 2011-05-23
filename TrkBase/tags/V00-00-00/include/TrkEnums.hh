//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkEnums.hh,v 1.10 2003/06/29 13:17:20 raven Exp $
//
// Description:
//     Dummy class (== poor man's namespace) to hold tracking enums.  Should 
//  move TrkDirection in here as well, but it's not worth the pain.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------

#ifndef TRKENUMS_HH
#define TRKENUMS_HH


namespace TrkEnums {

  enum TrkViewInfo {noView=-1,xyView=0, zView, bothView };
// flag possible instance meanings.  Note the range is _extremely_ restricted,
// so don't add values lightly.
  enum PackFlag{KalFit=0,KalConstraint=1,maxval=3};

};

#endif
