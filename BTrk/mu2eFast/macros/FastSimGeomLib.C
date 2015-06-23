#include <TString.h>
#include <Riostream.h>

#include <iostream>
#include <vector>

using namespace std;

TString MeasFileName="";
TString GeomFileName="";

void writeHeader(ofstream *outfile){

  (*outfile)<<"<?xml version=\"1.0\" encoding=\"UTF-8\" ?>"<<endl;
  (*outfile)<<"<edml>"<<endl;
  (*outfile)<<"    <included>"<<endl;

}

void writeFooter(ofstream *outfile){

  (*outfile)<<"    </included>"<<endl;
  (*outfile)<<"</edml>"<<endl;

}

ofstream * OpenMeasFile(TString fname){
  fname+="_Measures.xml";
  MeasFileName=fname;
  ofstream *measfile = new ofstream(fname.Data());
  writeHeader(measfile);
  (*measfile)<<"        <measures>"<<endl;
  (*measfile)<<endl;
  return measfile;
}

void CloseMeasFile(ofstream *measfile){
  (*measfile)<<endl;
  (*measfile)<<"        </measures>"<<endl;
  writeFooter(measfile);
  measfile->close();
  delete measfile;
}

ofstream * OpenGeomFile(TString fname){
  fname+="_Geom.xml";
  GeomFileName=fname;
  ofstream *geomfile = new ofstream(fname.Data());
  writeHeader(geomfile);
  (*geomfile)<<"        <detector name=\"Mu2e\">"<<endl;
  (*geomfile)<<"            <volume name=\"Tracker\">"<<endl;
  (*geomfile)<<endl;
  return geomfile;
}

void CloseGeomFile(ofstream *geomfile){
  (*geomfile)<<endl;
//   (*geomfile)<<"            </volume>"<<endl;
//   (*geomfile)<<"        </detector>"<<endl;
  writeFooter(geomfile);
  geomfile->close();
  delete geomfile;
}

void writeDetectorFile(TString fname){
  fname+=".xml";
  if (MeasFileName.IsNull() || GeomFileName.IsNull()) {
    cerr<<"You must open the Measure and Geometry files before!!!!"<<endl;
  }
  ofstream *detfile = new ofstream(fname.Data());
  writeHeader(detfile);
  (*detfile)<<endl;
  (*detfile)<<"        <include file=\""<<MeasFileName<<"\"  />"<<endl;
  (*detfile)<<"        <include file=\""<<GeomFileName<<"\"  />"<<endl;
  (*detfile)<<endl;
  writeFooter(detfile);
  detfile->close();
  delete detfile;
}

void writeDCHmeasure(ofstream *measfile, TString name, double cell_size=0.8073, double angle=0.0, double sTW=0.15e-6, double rms_par0=0.01,
 double rms_par1=0.0, double rms_par2=0.0, double rms_par3=0, double rms_par4=0.0, double rms_par5=0.0, double eff_par0=0.99,
 double eff_par1=1.00, double trunc_frac=0.8, double dedx_par1=0.0013202, double dedx_par2=1.0,
 double dedx_par3=-0.5){

  (*measfile)<<"            <device name=\""<<name<<"\""<<endl;
  (*measfile)<<"                type=\"DriftChamber\""<<endl;
  (*measfile)<<"                sensitiveTimeWindow=\""<<scientific<<sTW<<"\""<<endl;
  (*measfile)<<resetiosflags ( ios::scientific );
  (*measfile)<<"                rms_par0=\""<<rms_par0<<"\""<<endl;
  (*measfile)<<"                rms_par1=\""<<rms_par1<<"\""<<endl;
  (*measfile)<<"                rms_par2=\""<<rms_par2<<"\""<<endl;
  (*measfile)<<"                rms_par3=\""<<rms_par3<<"\""<<endl;
  (*measfile)<<"                rms_par4=\""<<rms_par4<<"\""<<endl;
  (*measfile)<<"                rms_par5=\""<<rms_par5<<"\""<<endl;
  (*measfile)<<"                eff_par0=\""<<eff_par0<<"\""<<endl;
  (*measfile)<<"                eff_par1=\""<<eff_par1<<"\""<<endl;
  (*measfile)<<"                cell_size=\""<<cell_size<<"\""<<endl;
  (*measfile)<<"                trunc_frac=\""<<trunc_frac<<"\""<<endl;
  (*measfile)<<"                dedx_par1=\""<<dedx_par1<<"\""<<endl;
  (*measfile)<<"                dedx_par2=\""<<dedx_par2<<"\""<<endl;
  (*measfile)<<"                dedx_par3=\""<<dedx_par3<<"\""<<endl;
  (*measfile)<<"                angle=\""<<angle<<"\"  />"<<endl;

}

void writeCyl(ofstream *outfile, TString name, int id, double zmin, double zmax, double radius, double thick,
 TString mat="", TString meas="", double gap=-1.0, double overlap=-1.0){

  (*outfile).setf(ios::fixed,ios::floatfield);
  (*outfile).precision(6);
  (*outfile)<<"                <cyl name=\""<<name<<"\" id=\""<<id<<"\" zmin=\""<<zmin<<"\" zmax=\""<<zmax;
  (*outfile)<<"\" radius=\""<<radius<<"\" thick=\""<<thick<<"\"";
  if (!mat.IsNull()) (*outfile)<<" mat=\""<<mat<<"\"";
  if (!meas.IsNull()) (*outfile)<<" meas=\""<<meas<<"\"";
  if (gap>0.0) (*outfile)<<" gap=\""<<gap<<"\"";
  if (overlap>0.0) (*outfile)<<" overlap=\""<<overlap<<"\"";
  (*outfile)<<"  />"<<endl;

}

void writeRing(ofstream *outfile, TString name, int id, double z, double lowradius, double hiradius, double thick,
 TString mat="", TString meas="", double gap=-1.0, double overlap=-1.0){

  (*outfile).setf(ios::fixed,ios::floatfield);
  (*outfile).precision(6);
  (*outfile)<<"                <ring name=\""<<name<<"\" id=\""<<id<<"\" z=\""<<z<<"\" lowradius=\""<<lowradius;
  (*outfile)<<"\" hiradius=\""<<hiradius<<"\" thick=\""<<thick<<"\"";
  if (!mat.IsNull()) (*outfile)<<" mat=\""<<mat<<"\"";
  if (!meas.IsNull()) (*outfile)<<" meas=\""<<meas<<"\"";
  if (gap>0.0) (*outfile)<<" gap=\""<<gap<<"\"";
  if (overlap>0.0) (*outfile)<<" overlap=\""<<overlap<<"\"";
  (*outfile)<<"  />"<<endl;

}

void writeVoxelization(ofstream *outfile, std::vector<double> &rbou, std::vector<int> &nphi, std::vector<double> &zbou){

  (*outfile)<<endl;
  (*outfile)<<"            </volume>"<<endl;
  (*outfile)<<"        </detector>"<<endl;

  (*outfile).setf(ios::fixed,ios::floatfield);
  (*outfile).precision(6);
  (*outfile)<<"<!-- now define the voxelization -->"<<endl;
  (*outfile)<<"        <config>"<<endl;
  (*outfile)<<"            <sect name=\"Tracker\">"<<endl;
  (*outfile)<<"                <param name=\"rbounds\"    type=\"vector\"        >"<<endl;
  (*outfile)<<"                   ";
  for (unsigned i=0; i<rbou.size(); i++){
    (*outfile)<<" "<<rbou.at(i);
  }
  (*outfile)<<endl;
  (*outfile)<<"                </param>"<<endl;
  (*outfile)<<"                <param name=\"nphivoxels\"    type=\"vector\"        >"<<endl;
  (*outfile)<<"                   ";
  for (unsigned i=0; i<nphi.size(); i++){
    (*outfile)<<" "<<nphi.at(i);
  }
  (*outfile)<<endl;
  (*outfile)<<"                </param>"<<endl;
  (*outfile)<<"                <param name=\"zbounds\"    type=\"vector\"        >"<<endl;
  (*outfile)<<"                   ";
  for (unsigned i=0; i<zbou.size(); i++){
    (*outfile)<<" "<<zbou.at(i);
  }
  (*outfile)<<endl;
  (*outfile)<<"                </param>"<<endl;
  (*outfile)<<endl;
  (*outfile)<<"            </sect>"<<endl;
  (*outfile)<<"        </config>"<<endl;

}
