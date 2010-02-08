// File and Version Information:
//      $Id: KalCodes.hh,v 1.7 2005/09/08 19:07:43 brownd Exp $
//
// Description:
//      define codes used in KalRep fitting
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 3/15/97
//------------------------------------------------------------------------
#ifndef KALCODES_HH
#define KALCODES_HH
// use a class as a namespace
class KalCodes {
public:
  enum kalFitFailCodes { dof=11,stops=12,diverge=13,matrix=14,processing=15,
			 notready=16,stubmatch=17,extendable=18,momentum=19,
			 endextend=20,inconsistent=21,nostop=22,cannotfullyextend=23,
			 badtraj=24};
  enum kalFitSuccessCodes {current=11,valid=12,alreadyextended=13};
  enum kalMakerFailCodes { makerep=30,fitresult=31,notkalrep=32,norep=33};
  enum kalMakerSuccessCodes { hypoexists=30};
  enum KalMiniSuccessCodes {unchanged=40};
  enum KalMiniFailCodes {nohots=40,nocache=41,mustbeactive=42};
};
#endif
