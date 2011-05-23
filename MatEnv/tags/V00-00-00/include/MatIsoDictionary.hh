//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatIsoDictionary.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatIsoDictionary (Isotope Dictionary)
//      Header file
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby  (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   October 15, 1998 - Talby : created
//-----------------------------------------------------------------------------

#ifndef MATISODICTIONARY_HH
#define MATISODICTIONARY_HH

//----------------------
// Base Class Headers --
//----------------------
#include <string>
#include <map>
#include "BaBar/BbrCollectionUtils.hh"

#include "MatEnv/MatIsotopeObj.hh"
#include "MatEnv/MatIsotopeList.hh"
//-------------------------------
// Collaborating Class Headers --
//-------------------------------

class MatIsoDictionary : public std::map<std::string*, MatIsotopeObj*, babar::Collection::PtrLess>
{

public:

// Constructor 
  MatIsoDictionary();
 
// Destructor
  virtual ~MatIsoDictionary();

  void FillIsoDict(MatIsotopeList* isolist);
};

#endif /* MATISODICTIONARY_HH */

