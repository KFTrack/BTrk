#ifndef BaBaR_Pdt_HH
#define BaBaR_Pdt_HH
//
// Wrapper around mu2e::ParticleDataTable class to make it look and feel
// like the BaBaR Pdt class.
//
// $Id:$
// $Author:$
// $Date:$
//

#include "BaBar/PdtPid.hh"
#include "BaBar/PdtEntry.hh"

//#include "ConditionsService/inc/ConditionsHandle.hh"
//#include "ConditionsService/inc/ParticleDataTable.hh"

class Pdt{
public:
  static const PdtEntry* lookup(PdtPid::PidType){
    return 0;
  }
  static double mass(PdtPid::PidType){
    return 0.13956;
  }
};

#endif
