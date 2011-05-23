//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: MatDBInfo.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//	Class MatDBInfo.  Implementation of MaterialInfo interface
//      using the database.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Dave Brown                      LBL
//
// Copyright Information:
//	Copyright (C) 1999		Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------

#ifndef MATDBINFO_HH
#define MATDBINFO_HH

#include "MatEnv/MaterialInfo.hh"

class DetMaterial;
class RecoMatFactory;
#include <string>
#include <map>
#include "BaBar/BbrCollectionUtils.hh"
using babar::Collection::PtrLess;

class MatBuildEnv;

class MatDBInfo : public MaterialInfo {
public:
  MatDBInfo();
  virtual ~MatDBInfo();
//  Find the material, given the name
  virtual const DetMaterial* findDetMaterial( const std::string& matName ) const;
// utility functions
private:
  DetMaterial* createMaterial( const std::string& dbName, 
			       const std::string& detMatName ) const;
  void declareMaterial( const std::string& dbName, 
			const std::string& detMatName );
// Cache of RecoMatFactory pointer
  RecoMatFactory* _genMatFactory;
// Cache of list of materials for DetectorModel
  std::map< std::string*, DetMaterial*, PtrLess > _matList;
// Map for reco- and DB material names
  std::map< std::string, std::string > _matNameMap; 
// function to cast-off const
  MatDBInfo* that() const {
    return const_cast<MatDBInfo*>(this);
  }
// allow MatBuildEnv to mess with me
  friend class MatBuildEnv;
  friend class MatBuildCoreEnv;
};

#endif
