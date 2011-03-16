//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatMaterialList.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatMaterialList (transient version)
//      Source file
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby  (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   June 18, 1998 - Talby : created
//-----------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
//----------------------
// C++ Headers --
//----------------------
#include <fstream>
#include <stdlib.h>
#include <assert.h>
#include <algorithm>
#include "BbrStdUtils/BbrCollectionUtils.hh"
using babar::Collection::DeleteObject;

//----------------------
// Base Class Headers --
//----------------------
#include "MatEnv/MatMaterialList.hh"
#include "ErrLogger/ErrLog.hh"
using std::fstream;
using std::ifstream;

//-------------------------------
// Collaborating Class Headers --
//-------------------------------


// Constructor to create an Material

MatMaterialList::MatMaterialList()
  : _vector(0)
{
}

MatMaterialList::MatMaterialList(const std::vector<MatMaterialObj*>& vector)
  : _vector(vector)
{
}

MatMaterialList::MatMaterialList(const std::string& materialsFile) 

{

// open input file materialsFile to read materials one by one 
  ifstream materials( materialsFile.c_str() );
  assert( materials.good() );
  if (materials.eof()) {
    ErrMsg(fatal) << "MatEnv/MatMaterialsList.data file empty!" << endmsg; 
  }

  std::string tagname;
  std::string fline;

//  Read, skipping comments
  tagname = "Materials_list";
  bool tag = false;
  while(!tag){
    do {  
      getline(materials, fline);
    } while (fline == "" && !materials.eof());
    tag = ( fline.find(tagname) != std::string::npos );
  }
  assert(tag);

// read Materials data (Name, Zeff, Aeff, nComp ...)
  std::string name;
  std::string cpname;
  std::string state;
  std::vector<double> Compweight;
  std::vector<std::string> Compname;
  std::vector<int> Compflag;
  int matentry = 0;
  int nbrcomp = 0;
  int iflag = 0;
  double density = 0.;
  double zeff = 0.;
  double aeff = 0.;
  double weight = 0.;
  double radlen = 0.;
  double intlen = 0.;
  double refindex = 0.;
  double temperature = 0.;
  double pressure = 0.;
  materials >> name;
  while( !materials.eof())
    {
      if(name.find("#") != std::string::npos) {
	do {  
	  getline(materials, fline);
	} while (fline == "" && !materials.eof());
	materials >> name;
      } else {
	materials >> density >> zeff >> aeff >> nbrcomp;
	if(nbrcomp != 0) {
	  for (int i=0; i<abs(nbrcomp); i++) 
	    { 
	      materials >> weight >> cpname >> iflag;
	      Compflag.push_back(iflag);
	      Compweight.push_back(weight);
	      Compname.push_back(cpname);
	    }

	} else {
          Compflag.push_back(0);
          Compweight.push_back(0);
          Compname.push_back(" ");
	}            
	materials >> radlen >> intlen >> refindex >> temperature 
		  >> pressure >> state;

	MatMaterialObj* matObj = new MatMaterialObj();
	matObj->setName(name);
	matObj->setDensity(density);
	matObj->setZeff(zeff);
	matObj->setAeff(aeff);
	matObj->setNbrComp(nbrcomp);
	matentry = Compflag.size();
	for (size_t idx=matentry-abs(nbrcomp); idx<matentry; idx++)
	  {
	    matObj->setIflg(Compflag[idx]);
	    matObj->setWeight(Compweight[idx]);
	    matObj->setCompName(Compname[idx]); 
	  }
	matObj->setRadLength(radlen);
	matObj->setIntLength(intlen);
	matObj->setRefIndex(refindex);
	matObj->setTemperature(temperature);
	matObj->setPressure(pressure);
	matObj->setState(state);

	_vector.push_back(matObj);

//      _vector.push_back(new MatMaterialObj(name, density, zeff, 
//                         aeff, nbrcomp, Compflag, Compweight, 
//                         Compname, radlen, intlen, refindex, 
//                         temperature, pressure, state));
      materials >> name;
      }
    }  

}

MatMaterialList::~MatMaterialList() 
{
  std::for_each(_vector.begin(), _vector.end(), DeleteObject());
  _vector.clear();   
}






