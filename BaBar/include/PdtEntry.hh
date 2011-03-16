#ifndef BaBaR_PdtEntry_HH
#define BaBaR_PdtEntry_HH
//
// Wrapper around mu2e::ParticleDataTable class to make it look and feel
// like the BaBaR PdtEntry class.
//
// $Id:$
// $Author:$
// $Date:$
//

#include <string>

class PdtEntry{
public:
  std::string name()const {return "";};
  double mass() const{return 0.;};
};

#endif
