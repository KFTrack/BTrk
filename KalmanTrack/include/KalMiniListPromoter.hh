//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: KalMiniListPromoter.hh,v 1.1 2005/04/04 19:54:54 brownd Exp $
//
// Description:
//	Module KalMiniListPromoter.  Promote the tracks on the
//      specified list to 'hots' mode.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//     David Brown  April 1, 2005
//
// Copyright Information:
//	Copyright (C) 2005		Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------

#ifndef KalMiniListPromoter_HH
#define KalMiniListPromoter_HH

//----------------------
// Base Class Headers --
//----------------------
#include "Framework/AppModule.hh"
#include "AbsParm/AbsParmIfdStrKey.hh"
#include "AbsParm/AbsParmNamedValue.hh"

class KalMiniListPromoter : public AppModule {

//--------------------
// Instance Members --
//--------------------

public:
// specify which fits to promote.
  enum fitHypos { defhypo=0,allexisting};
// Constructors
  KalMiniListPromoter( const char* const theName, const char* const theDescription );
  KalMiniListPromoter( const KalMiniListPromoter& other, const char* const theName );
// clone function (with covariant return)
  virtual KalMiniListPromoter* clone(const char* cloneName);
// Destructor
  virtual ~KalMiniListPromoter( );
// Operations
  virtual AppResult beginJob( AbsEvent* anEvent );
  virtual AppResult event( AbsEvent* anEvent );
private:
  AbsParmIfdStrKey _trklist; // track list
  AbsParmNamedValue<fitHypos> _hypos; // hypos to promote
};
#endif
