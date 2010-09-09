//
// BaBar includes
#include "BaBar/BaBar.hh"
#include "Framework/AppFileName.hh"
//specific to this code, taken from mu2e framework
#include "mu2eFast/mu2eDSField.hh"
#include "mu2eFast/BFMap.hh"
#include "mu2eFast/DiskRecord.hh"
#include "mu2eFast/Container3D.hh"
#include "mu2eFast/MinMax.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

// Includes from C ( needed for block IO ).
#include <fcntl.h>
#include <sys/stat.h>

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <limits>
#include <cmath>
#include <assert.h>

using namespace std;

// DS orgin in mu2e official coordinates (in mm!)
Hep3Vector
mu2eDSField::trackerCenterInMu2eCoordinates = Hep3Vector(-3904.,0.,10200.);

// Fastsim origin
Hep3Vector
mu2eDSField::trackerCenterInFastSimCoordinates = Hep3Vector(0.0,0.0,0.0);

mu2eDSField::mu2eDSField(const std::string& fname,double dfactor) : _dfactor(dfactor){
  // parameters for the DS field map taken from geom_01.txt
  int nx = 50;
  int ny = 25;
  int nz = 438;
  bool warnIfOutside = false;
  std::string filename = AppFileName(fname.c_str()).pathname();
  std::string key("DS");
  // Create an empty map
  Hep3Vector maporigin;
  _fieldmap = mu2e::BFMap(key,maporigin,nx,ny,nz,warnIfOutside);
  
  // Fill the map.
  readGMCMap( filename, _fieldmap );
  
  // set nominal field
  _bnom = _fieldmap.getBField(trackerCenterInMu2eCoordinates).z();
  // decide if we need to scale distortions
  _distort = _dfactor != 1.0;
}

mu2eDSField::~mu2eDSField(){}

// BaBar interface.  Note we have to change units here to the BaBar conventions
Hep3Vector
mu2eDSField::bFieldVect (const HepPoint &point)const {
  static Hep3Vector nomfield(0.0,0.0,_bnom);
// convert units (cm to mm) and origins to meco conventions
  Hep3Vector mecopoint(10*(point.x()-trackerCenterInFastSimCoordinates.x()) + trackerCenterInMu2eCoordinates.x(),
  10*(point.y()-trackerCenterInFastSimCoordinates.y()) + trackerCenterInMu2eCoordinates.y(),
  10*(point.z()-trackerCenterInFastSimCoordinates.z()) + trackerCenterInMu2eCoordinates.z());
  Hep3Vector mfield = _fieldmap.getBField(mecopoint);
  if(_distort){
// subtract out the nominal field, scale the difference, then add back the nominal
    mfield = nomfield + _dfactor*(mfield-nomfield);
  }
  return mfield;
}
  
double
mu2eDSField::bFieldNominal()const {
  return _bnom;
}


//
// Read one magnetic field map file in MECO GMC format.
//
// This does a 2 pass operation"
// 1) Pass 1:
//      Read the input file into a temporary image in memory.
//      Find the min and max values of the grid points.
//      A the end of this pass, compute the grid spacing.
// 2) Pass 2:
//      Fill the 3D array from the image in memory.
//
void
mu2eDSField::readGMCMap( const string& filename, mu2e::BFMap& bfmap ){

// Open the input file.
  int fd                                    = open( filename.c_str(), O_RDONLY );
  if ( !fd ) {
    cerr 
      << "Could not open file containing the magnetic filed map for: "
      << bfmap.getKey() << "\n"
      << "Filename: " 
      << filename
      << "\n" << endl;
  }

     // Compute number of records in the input file.
  const int nrecords                        = computeArraySize(fd,filename);

     // Image of the file in memory.
  vector<DiskRecord> data(nrecords, DiskRecord());

     // Read file into memory.
  const int nbytes                          = nrecords*sizeof(DiskRecord);
  ssize_t s                                 = read( fd, &data[0], nbytes );
  if ( s != nbytes ) {
    if ( s == -1 ){
      cerr
        << "Error reading magnetic field map: " 
        << bfmap.getKey() << "\n"
        << "Filename: " 
        << filename
        << "\n";
    } else{
      cerr
        << "Wrong number of bytes read from magnetic field map: " 
        << bfmap.getKey() << "\n"
        << "Filename: " 
        << filename
        << "\n";
    }
  }

     // Tool to find min and max values of grid points.
  MinMax mmX, mmY, mmZ;

     // Offset needed to put this map into the Mu2e coordinate system.
     // ( Origin at center of TS ).
  const CLHEP::Hep3Vector& offset           = bfmap.origin();

     // Collect distinct values of (X,Y,Z) on the grid points.
  set<float> X, Y, Z;

     // Multiply by this factor to convert from kilogauss to tesla.
  double ratio                              = CLHEP::kilogauss/CLHEP::tesla;

     // For the image in memory:
     //   1) Transform into the correct set of units.
     //   2) Find min/max of each dimension.
     //   3) Collect unique values of (X,Y,Z) of the grid points.
  for ( vector<DiskRecord>::iterator i      = data.begin();
  i != data.end(); ++i ){

       // Modify in place.
    DiskRecord& r                           = *i;

       // Unit conversion: from (cm, kG) to (mm,T).
    r.x  *= CLHEP::cm; 
    r.y  *= CLHEP::cm; 
    r.z  *= CLHEP::cm;
    r.bx *= ratio; 
    r.by *= ratio; 
    r.bz *= ratio; 

       // Re-centering. - not sure if we really want this here?
    r.x += offset.x(); 
    r.y += offset.y(); 
    r.z += offset.z();

       // The one check I can do.
    if ( r.head != r.tail ){
      cerr
        << "Error reading magnetic field map.  "
        << "Mismatched head and tail byte counts at record: " << data.size() << "\n"
        << "Could not open file containing the magnetic filed map for: "
        << bfmap.getKey() << "\n"
        << "Filename: " 
        << filename
        << "\n";
    }

       // Update min/max information.
    mmX.compare(r.x);
    mmY.compare(r.y);
    mmZ.compare(r.z);

       // Populate the set of all unique grid values.
    X.insert(r.x);
    Y.insert(r.y);
    Z.insert(r.z);
  }

     // Expected grid dimentsions.
  const int nx                              = bfmap.nx();
  const int ny                              = bfmap.ny();
  const int nz                              = bfmap.nz();

     // Cross-check that the grid read from the file has the size we expected.
     // This is not really a perfect check since there could be round off error
     // in the computation of the grid points.  But the MECO GMC files were written 
     // in a way that this test works.
  if ( X.size() != nx ||
    Y.size() != ny ||
  Z.size() != nz     ){
    cerr
      << "Mismatch in expected and observed number of grid points for BField map: " 
      << bfmap.getKey() << "\n"
      << "From file: " 
      << filename
      << "\n"
      << "Expected/Observed x: " << nx << " " << X.size() << "\n"
      << "Expected/Observed y: " << ny << " " << Y.size() << "\n"
      << "Expected/Observed z: " << nz << " " << Z.size() << "\n";
  }

     // Cross-check that we did not find more data than we have room for.
  if ( data.size() > nx*ny*nz-1){
    cerr
      << "Too many values read into the field map: " 
      << bfmap.getKey() << "\n"
      << "From file: " 
      << filename
      << "\n"
      << "Expected/Observed size: " << nx*ny*nz << " " << data.size() << "\n";
  }

     // Tell the magnetic field map what its limits are.
  bfmap.setLimits( mmX.min(), mmX.max(),
    mmY.min(), mmY.max(),
    mmZ.min(), mmZ.max() );

     // Store grid points and field values into 3D arrays
  for (vector<DiskRecord>::const_iterator i = data.begin(), e=data.end();
  i != e; ++i){

    DiskRecord const& r(*i);

       // Find indices corresponding to this grid point.
       // By construction the indices must be in bounds ( we set the limits above ).
    std::size_t ix                          = bfmap.iX(r.x);
    std::size_t iy                          = bfmap.iY(r.y);
    std::size_t iz                          = bfmap.iZ(r.z);

       // Store the information into the 3d arrays.
    bfmap.grid().set (ix, iy, iz, CLHEP::Hep3Vector(r.x,r.y,r.z));
    bfmap.field().set(ix, iy, iz, CLHEP::Hep3Vector(r.bx,r.by,r.bz));
    bfmap.setDefined(ix, iy, iz, true);
    
  }

  return;

}

// Compute the size of the array needed to hold the raw data of the field map.
int
mu2eDSField::computeArraySize( int fd, const string& filename ){

  // Get the file size, in bytes, ( info.st_size ).
  struct stat info;
  fstat( fd, &info);

  // Check that an integral number of records fits in the file.
  int remainder = info.st_size % sizeof(DiskRecord);
  if ( remainder != 0 ){
    cout 
      << "Field map file does not hold an integral number of records: \n" 
      << "Filename:  " << filename << "\n"
      << "Size:      " << info.st_size << "\n"
      << "Remainder: " << remainder << "\n";
    assert(2==1);
  }

  // Compute number of records.
  int nrecords  = info.st_size / sizeof(DiskRecord);
  return nrecords;
}
