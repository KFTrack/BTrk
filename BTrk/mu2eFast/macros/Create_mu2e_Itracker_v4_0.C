#define firstsenserad      38.1972000
#define firstlayercellwd    0.6000000
#define ncell             400
#define swirediam           0.0020000
#define nbasefwirepercell   4
#define fwirediam           0.0040000
#define nlayer             38
#define averageangle        0.1500000

#define averageresol        0.0100000     //cm
#define vdrift              3.0000000e+6  //cm/s

#define legth             100.0000000
#define zoffset           -50.0
	    
#define innerwallthickness  0.0266000     // equivalent in grafhite of the inner wall sandwitch
#define endplatethickness   0.0445000     // equivalent in Aluuminium
#define outerwallthickness  1.0000000

#define voxelsafetydim      0.01

const TString gasMix        ="IT-gas1";
const TString fwMaterial    ="IT-Fwire";
const TString swMaterial    ="IT-Swire";
const TString EPeq_Material ="IT-Aluminum";
const TString IWeq_Material ="IT-Graphite";
const TString OWeq_Material ="CFiber";


void Create_mu2e_Itracker_v4_0(TString fileSuffix = "mu2e_Itracker"){

  FileStat_t  finfoSRC, finfoLIB;
  if (gSystem->GetPathInfo("mu2eFast/macros/FastSimGeomLib_C.so",finfoLIB)!=0){
    gROOT->ProcessLine(".L mu2eFast/macros/FastSimGeomLib.C++");
  }
  else {
    gSystem->GetPathInfo("mu2eFast/macros/FastSimGeomLib.C",finfoSRC);
    if (finfoSRC.fMtime>finfoLIB.fMtime) gROOT->ProcessLine(".L mu2eFast/macros/FastSimGeomLib.C++");
    else gSystem->Load("mu2eFast/macros/FastSimGeomLib_C.so");
  }  

  
  ofstream *measf = OpenMeasFile(fileSuffix);
  ofstream *geomf = OpenGeomFile(fileSuffix);

//----- building the geometry -----
  double rScaleCoeff    = (firstsenserad+0.5*firstlayercellwd)/(firstsenserad-0.5*firstlayercellwd);
  double gaszmin        = zoffset;
  double gaszmax        = legth + zoffset;
  double fwirethickness = TMath::Pi()*pow(0.5*fwirediam,2)/fwirediam;
  double swirethickness = TMath::Pi()*pow(0.5*swirediam,2)/swirediam;
  double fweqlength     = fwirediam*nbasefwirepercell;
  double cellwd         = firstlayercellwd;
  double swrad          = firstsenserad;
  
  std::vector<double>   voxRbound;
  std::vector<int>      voxNphi;
  std::vector<double>   voxZbound;
 
  int nomeasid;
  int idGasDown  = 1+nlayer;
  int idGasUp    = 1+2*nlayer;
  int idCenterFw = 1+3*nlayer;
  int idSw       = 1+4*nlayer;
  int idOuterFw  = 1+5*nlayer;
  
  nomeasid = 1+6*nlayer;
  
  double innerwallrad   = firstsenserad-0.5*(/*fwirediam+*/firstlayercellwd+innerwallthickness);
  writeCyl(geomf,"dch-InnerCyl",++nomeasid,gaszmin,gaszmax,innerwallrad,
    innerwallthickness,IWeq_Material);   //Inner Wall

  voxRbound.push_back(innerwallrad-0.5*innerwallthickness);
  voxNphi.push_back(1);
  voxRbound.push_back(innerwallrad+0.5*innerwallthickness);

  //nomeasid=nlayer;

  writeCyl(geomf,"dch-inner-field-wires",0/*++nomeasid*/,gaszmin,gaszmax,(firstsenserad-0.5*(firstlayercellwd-fwirediam)),
    fwirethickness,fwMaterial, "", (1.0-fweqlength/cellwd));   //innermost fw

  TString measureName="";
  double sangle=0.0;
  
  double alfCellwd;
  
  for (int il=1; il<=nlayer; il++){
     measureName="Stereo";
     if (il%2 ==0) {
       measureName+=Form("%i+",il);
       sangle=averageangle;
     }
     else {
       measureName+=Form("%i-",il);
       sangle=-1.0*averageangle;
     }
     
     //cout<<"drop "<<TMath::Sqrt(swrad*swrad + 0.250*pow(legth*TMath::Tan(0.150),2)) - swrad<<endl;
     
     alfCellwd = 0.5*cellwd;
     
     idGasDown  ++;
     idGasUp    ++;
     idCenterFw ++;
     idSw       ++;  
     idOuterFw  ++;

     writeDCHmeasure(measf,measureName,cellwd,sangle,alfCellwd/vdrift,averageresol);
     (*measf)<<endl;

     writeCyl(geomf,"dch-Gas_down",idGasDown/*il*/,gaszmin,gaszmax,swrad-0.5*alfCellwd,alfCellwd,gasMix);  //gas of the alf cell before the sense  wires
     writeCyl(geomf,"dch-meas",il,gaszmin,gaszmax,swrad,cellwd,"",measureName);  //measure layer
     writeCyl(geomf,"dch-Gas_up",idGasUp/*il*/,gaszmin,gaszmax,swrad+0.5*alfCellwd,alfCellwd,gasMix);  //gas of the alf cell after the sense  wires
     writeCyl(geomf,"dch-center-field-wires",idCenterFw/*++nomeasid*/,gaszmin,gaszmax,swrad-0.01,fwirethickness,fwMaterial, "",
       (1.0-fwirediam/cellwd));  //central fw
     writeCyl(geomf,"dch-sense-wires",idSw/*++nomeasid*/,gaszmin,gaszmax,swrad+0.01,swirethickness,swMaterial, "",
       (1.0-swirediam/cellwd));  //sw
       writeCyl(geomf,"dch-out-field-wires",idOuterFw/*++nomeasid*/,gaszmin,gaszmax,swrad+alfCellwd-0.5*fwirediam,fwirethickness,fwMaterial, "",
       (1.0-fweqlength/cellwd));  //upper fw

     voxNphi.push_back(1);  //ncell
     voxRbound.push_back(swrad+cellwd);
     swrad*=rScaleCoeff;
     cellwd*=rScaleCoeff;
  }

  
  //nomeasid = nlayer*5+1;
  
  double outerwallrad    = (swrad+cellwd)/rScaleCoeff + 0.5*outerwallthickness;
  writeCyl(geomf,"dch-OuterCyl",++nomeasid,gaszmin,gaszmax,outerwallrad,outerwallthickness,OWeq_Material);
  innerwallrad-=0.5*innerwallthickness;
  outerwallrad+=0.5*outerwallthickness;
  
  voxNphi.push_back(1);
  voxRbound.push_back(outerwallrad+voxelsafetydim);

  writeRing(geomf,"dch-Endplate",++nomeasid,gaszmin-endplatethickness-voxelsafetydim,innerwallrad,outerwallrad,
    endplatethickness,EPeq_Material);
  writeRing(geomf,"dch-Endplate",++nomeasid,gaszmax+endplatethickness+voxelsafetydim,innerwallrad,outerwallrad,
    endplatethickness,EPeq_Material);

  voxZbound.push_back(gaszmin-endplatethickness-endplatethickness-voxelsafetydim);
  voxZbound.push_back(gaszmin-voxelsafetydim);
  voxZbound.push_back(gaszmax+voxelsafetydim);
//   double zstep = 0.1000000;
//   double tempz;
//   tempz = gaszmin;
//   while(true){
//     voxZbound.push_back(tempz);
//     tempz+=zstep;
//     if (tempz>=gaszmax) break;
//   }
  voxZbound.push_back(gaszmax+endplatethickness+endplatethickness+voxelsafetydim);

//   for (int i=0; i<(voxRbound.size()-1); i++)
//     voxNphi.push_back(1);

  writeVoxelization(geomf,voxRbound,voxNphi,voxZbound);

//----- 

  CloseMeasFile(measf);
  CloseGeomFile(geomf);
  writeDetectorFile(fileSuffix);

}

