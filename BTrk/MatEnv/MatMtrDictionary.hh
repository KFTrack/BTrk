//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatMtrDictionary.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatMtrDictionary (Material Dictionary)
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

#ifndef MATMTRDICTIONARY_HH
#define MATMTRDICTIONARY_HH

//----------------------
// Base Class Headers --
//----------------------
#include <string>
#include <map>
#include "BTrk/BaBar/BbrCollectionUtils.hh"

#include "BTrk/MatEnv/MatMaterialObj.hh"
#include "BTrk/MatEnv/MatMaterialList.hh"
//-------------------------------
// Collaborating Class Headers --
//-------------------------------

class MatMtrDictionary : public std::map<std::string*, MatMaterialObj*, babar::Collection::PtrLess>
{

public:

// Constructor 
  MatMtrDictionary();

// Destructor
  virtual ~MatMtrDictionary();

  void FillMtrDict(MatMaterialList* mtrlist);

};

#endif /* MATMTRDICTIONARY_HH */




