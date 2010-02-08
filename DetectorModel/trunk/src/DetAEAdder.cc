// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetAEAdder.cc,v 1.1 2002/01/15 21:34:23 brownd Exp $
//
//  Description:  DetAEAdder
//
// Copyright Information:
//	Copyright (C) 2002	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 1/15/2002
//------------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "DetectorModel/DetAEAdder.hh"
#include "DetectorModel/DetAlignElem.hh"
#include "CLHEP/Matrix/SymMatrix.h"

bool
DetAEAdder::combine(const DetAlignElem& inone,
		    const DetAlignElem& intwo,
		    DetAlignElem& out) const
{
// user operator *= to set the parameter vector
  out = inone;
  out *= intwo;
// set the covariance
  out.setCovariance(intwo.parameterCovariance());
  return true;
}
