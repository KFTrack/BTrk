//
//  Code for the DetElem class.
//
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/ErrLog.hh"
#include <assert.h>
#include <vector>
#include <string>
#include "BTrk/BbrGeom/Transformation.h"
#include "BTrk/DetectorModel/DetElem.hh"
#include "BTrk/DetectorModel/DetAlignElem.hh"
#include "BTrk/DetectorModel/DetIntersection.hh"
#include "BTrk/DetectorModel/GnuPlot.hh"
#include "BTrk/DetectorModel/DetMaterial.hh"
using std::endl;
using std::ofstream;
using std::ostream;


//
//  default constructor
//
DetElem::DetElem()
: _ielem(0),_ename("Unknown"),  _dtype(0)
{;}
//
//  Real constructor
//
DetElem::DetElem(const DetType* type,const char* name,int id)
: _ielem(id),_ename(name),_dtype(type)
{;}
DetElem::~DetElem(){
}
//
//  print output;
//
void
DetElem::print(ostream& os) const {
  os << "Detector Element # " << _ielem << " " << _ename << endl;
}
//
void
DetElem::printAll(ostream& os) const {
  os << "Detector Element # " << _ielem << " " << _ename << " of ";
  _dtype->print(os);
}
//
//  Apply the global alignment correction.  This uses the
//  transform describing the local coordinate system
//  WRT the BaBar coordinate system.
//
void
DetElem::applyGlobal(const DetAlignElem& global){
  _etrans->transform(global.transform());
  updateCache();
}
//
//  Apply the local alignment correction.  This is the transform
//  describing the physical position of the wafer WRT its nominal position
//  in the nominal global coordinate system.
//
void
DetElem::applyLocal(const DetAlignElem& local){
  assert(local.elementNumber() == elementNumber());
  *(_etrans) *= local.transform();
  updateCache();
}
//
//  Remove the alignment; just invert the multiplication order and
//  the transforms themselves relative to the above.
//
void
DetElem::removeGlobal(const DetAlignElem& global){
  _etrans->transform(global.inverseTransform());
  updateCache();
}
void
DetElem::removeLocal(const DetAlignElem& local){
  assert(local.elementNumber() == elementNumber());
  *(_etrans) *= local.inverseTransform();
  updateCache();
}
//
//  dummy implementation of cache updating; used by the above
//  alignment functions (among other possibilities)
//
void
DetElem::updateCache(){
}
//
// Default implementation for the gnu plotting function
void
DetElem::gnuPlot( GnuPlot* gp ) const {
  ofstream* o  = gp->output;
  gp->stlPoints.clear();
  physicalOutline(gp->stlPoints);
  int ncorn = gp->stlPoints.size();
  if( ncorn > 0 ) {
    for(int icorn=0;icorn<ncorn;icorn++){
      *o << gp->stlPoints[icorn].x() << " "
	 << gp->stlPoints[icorn].y() << " "
	 << gp->stlPoints[icorn].z() << endl;
    }
    *o << endl;
  }

}
//
//  Implement dummy outline function.
//
void
DetElem::physicalOutline(std::vector<HepPoint>& pvec) const {
  pvec.clear();
}
//
//  Default implementation of material information function.  This should
//  be good enough for many subclasses.
//
void
DetElem::materialInfo(const DetIntersection& dinter,
		      double momentum,
		      TrkParticle const& tpart,
		      double& deflectRMS,
		      double& pFractionRMS,
		      double& pFraction,
		      trkDirection dedxdir) const {
  if(momentum>0.0){
//  Get the material
    const DetMaterial& mat = material(dinter);
//
//  The path length gives the effective thickness through the material
//
    double thickness = dinter.pathLength();
//
//  Pass everything to the material, to get the answers
//
    deflectRMS = mat.scatterAngleRMS(momentum,thickness,tpart);
    double echange = dedxdir == trkOut ? mat.energyLoss(momentum,thickness,tpart)
      : mat.energyGain(momentum,thickness,tpart);
    double elossRMS = mat.energyLossRMS(momentum,thickness,tpart);
//
//  Convert to (dimensionless) momentum fractions
//
    double energy = DetMaterial::particleEnergy(momentum,tpart);
    double newenergy = energy+echange;
    pFraction = DetMaterial::particleMomentum(newenergy,tpart)/momentum - 1.0;
    pFractionRMS = elossRMS*energy/(momentum*momentum);
  } else {
    deflectRMS = 1.0;
    pFraction = -1.0;
    pFractionRMS = 1.0;
  }
}
//
//  Default implementation of the material function, assuming that the
//  element is homogenous.  Non-homogenous element classes must
//  overwrite this.
//
const DetMaterial&
DetElem::material(const DetIntersection&) const {
  return detectorType()->material(0);
}

bool
DetElem::operator == (const DetElem& other) const {
  return _ielem == other._ielem && _ename == other._ename;
}

bool
DetElem::match(const DetAlignElem& ae) const {
  return elementNumber() == ae.elementNumber() &&
// align elem name might be truncated
    ( elementName().find(ae.elementName()) != std::string::npos );
}

bool 
DetElem::reIntersect(const Trajectory* traj,DetIntersection& dinter) const {
// update the intersection trajectory; global length stays the same
  dinter.delem = this;
  dinter.trajet = traj;
  return true;
}

