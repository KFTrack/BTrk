// -----------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetSurfaceRpcElem.hh,v 1.1 2008/05/14 16:23:50 fransham Exp $
//
//  Description:
//      A small class that inherits from DetSurfaceElem and is part of the Rpc
//
//  Authors: 
//	K. Fransham May 2008
//
//  Modifications:
//
//------------------------------------------------------------------------------
#include "BTrk/DetectorModel/DetSurfaceElem.hh"
class DetSurfaceRpcElem : public DetSurfaceElem {
public:
  DetSurfaceRpcElem(DetSurfaceType* type,const char* name,int id,const DetSurface& surf):
    DetSurfaceElem(type, name, id, surf) {}
};
