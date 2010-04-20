#ifndef PDTPID_HH
#define PDTPID_HH
#include "PDT/PdtLund.hh"

class PdtPid
{
public:
  enum {nPidType = 5};

  enum PidType
  {
    null = -1,
    electron = 0,
    muon = 1,
    pion = 2,
    kaon = 3,
    proton = 4
  };
  enum {nPidNeutralType = 5};

  enum PidNeutralType
  {
    none = -1,
    gamma = 0,
    pi0   = 1,
    K0L   = 2, 
    neutron = 3,
    anti_neutron = 4
  };
protected:
  
  friend class Pdt;
};

#endif
