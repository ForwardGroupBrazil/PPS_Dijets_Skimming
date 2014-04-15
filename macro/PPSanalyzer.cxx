#include "TROOT.h"
#include "TBox.h"
#include "TEllipse.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TDirectory.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TPaletteAxis.h"
#include "TPolyLine.h"
#include "TStyle.h"
#include "TText.h"
#include "TTree.h"
#include "TF1.h"
#include <iostream>
#include <string>
#include "PPSSpectrometer.h"
#include "TString.h"
#include "TLegend.h"


//Prototypes
void createCanvas();
void setBeamParameters(int pos, double x,double y,double rms_x,double rms_y);
void setDetectorParameters(float det_w = 22., float det_h=22.,float det_z1 = 203.8,float det_z2=215.,float det_p=15.,std::string br="Reco");
void drawHitmap(int norm=0);
void drawHitmap(std::string arm,int norm=0,float xi_min=0,float xi_max=0,float xmin=-30,float xmax=0,float ymin=-15, float ymax=15);
void draw_t_xi_acceptance(const std::string& arm, int norm=100, int nbins=100,double lot=0,double hit=4,double loxi=0.001,double hixi=0.2,std::string fmt="");
void draw_t_xi_acceptance(int norm = 100, int nbins=100,double lot=0,double hit=4.,double loxi=0.001,double hixi=0.2,string fmt="");
void draw_t_xi_rejection(int norm = 100, int nbins=100,double lot=0,double hit=4.,double loxi=0.001,double hixi=0.2,string fmt="");
void draw_totemhitmap(double RPxOffset=2.,double RPyOffset=4.,double xmin=-30.,double xmax=10.,double ymin=-20,double ymax=20);
void do_certification(double min_t=0,double max_t=7, double min_xi=0., double max_xi=0.3);
void draw_t_correlation(const char* ,double ,double );
void draw_xi_correlation(const char*,double ,double );
void draw_y_vs_mass_acceptance(double norm = 100., int nbins = 100, double mmin = 100., double mmax=500, double ymin=0, double ymax=2.,double ecms=13000.);
void draw_t_vs_mass_acceptance(const char*,double norm=100.,int nbins=100, double mmin=100, double mmax=1000,double tmin=0., double tmax=5.0,double ecms=13000.);
void draw_mass_acceptance(TCut="",const char* opt="", double ecms=13000.,double norm = 100, int nbins=100, double mmin=100, double mmax=1000);
void draw_centralmass(double ecms=13000, int nbins=100, double mmin=100, double mmax=1000);
void draw_mass_resolution(double mass,double ecms=13000., const char* br="Sim",const char* cut = "",int nbins=100, double d_mass=10);
void setZPosition(TTree* t);
void drawDetectorShape();
void drawBeamEnvelope(double sig);
void drawQuarticOccupancy(const string&, double xmin=-25, double xmax=0, double ymin=-12.5, double ymax=12.5);
void drawDiamondOccupancy(const string&, double xmin=-25, double xmax=0, double ymin=-12.5, double ymax=12.5);
void drawToFOccupancy(const string&, double xmin=-25, double xmax=0, double ymin=-12.5, double ymax=12.5);
void setToFParameters(std::string,double x);
void drawToFTotalMultiplicity(const std::string&);
void drawToFCellMultiplicity(const std::string&);
void drawToFHitPerCellMultiplicity(const std::string&);
void get_t_atMCframe(double cr_ang);

void setCollimator(const string& tcl, const float& pos);
void setDataFormat();
double ToFXFromCell(int celid);
double ToFYFromCell(int celid);
int  ToFCellId(double x, double y);
bool isValidCellId(int celid);
void setToFGeometry(std::string);
void setQuartic();
void setDiamond();
// Globals
namespace PPSanalyzer {
TCut detF;
TCut detB;
TCut det1F;
TCut det2F;
TCut det1B;
TCut det2B;
TCut ToFF;
TCut ToFB;
TCut momF;
TCut momB;
TCut TCL_F;
TCut TCL_B;
TCut stoppedF;
TCut stoppedB;
float detWidth;
float detHeight;
float detXPos1;  // first trakcer station X position (from beam in mm)
float detXPos2;  // second trakcer station X position (from beam in mm)
float detXPosSigma; // detectors distance from beam in sigmas
float detZPos1;
float detZPos2;
float ToFZPos;
bool fSingleArmAcceptance;
bool fToFAcceptance;
// fCheckPointIndex in the CheckPoint vector to look at for acceptances (GASTOF???)
int  fCheckPointIndex=0;
string fTarget="";
string fTreeName="";
string fDataFormat="";
string fRecoBranchName="";
string fSimBranchName="";
string fGenBranchName="";
bool   fInitialized=false;
bool   fRestrictedGen=false;
map<int,double> BeamX;
map<int,double> BeamY;
map<int,double> BeamX_RMS;
map<int,double> BeamY_RMS;
// ToF detector
std::string ToFGeometryType="";
double      ToFXPos;
vector<pair<double,double> > ToFCellColumn;
vector<pair<double,double> > ToFCellRow;
const double RPWindowThickness=0.3; 
const int   NYCell = 5; // number of cell in Y for Quartic
const int   NXCell = 5; // number of cell in X for Quartic
const int   NXCellLowerRow = 1; // number of cell in X for lower row for Diamond
const int   NXCellUpperRow = 1; // number of cell in X for upper row for Diamond
const int   NXCellCentralRow = 14; // number of cell in X for central row for Diamond
const int   DiamondNYCell = 3; // number of cell in Y for Diamond
const float celws[] = {0.3 , 0.4 , 0.45 , 0.55 , 0.65 , 0.75 , 0.9 , 1.0 , 1.3 , 1.7 , 2 , 2.4 , 3 , 4.6};
const float DiamondCellHeight = 5.;
const float QuarticCellHeight = 3.;
const float QuarticCellWidth = 3.;
const float DiamondUpperCellW = 10;
const vector<float> DiamondCentralCellW(celws,celws + sizeof(celws) / sizeof(float));
const float DiamondLowerCellW = 10;
};

using namespace PPSanalyzer;

void setDataFormat() {
     TTree *t = (TTree*)gROOT->FindObject("T");
     if (t) {
        fDataFormat = "PPS";
        fRecoBranchName="Reco";
        fSimBranchName="Sim";
        fGenBranchName="Gen";
        fTreeName="T";
        fInitialized=true;
        std::cout << "Selected PPS private format"<<endl;
     }
     else {
        t = (TTree*)gROOT->FindObject("Events");
        if (t) fDataFormat = "CMSSW";
        fRecoBranchName="PPSSpectrometer_ppssim_PPSReco_PPS.obj";
        fSimBranchName="PPSSpectrometer_ppssim_PPSSim_PPS.obj";
        fGenBranchName="PPSSpectrometer_ppssim_PPSGen_PPS.obj";
        fTreeName="Events";
        fInitialized=true;
        std::cout << "Selected CMSSW format"<<endl;
     }
     if (fDataFormat=="") {
        std::cout << "Environment not properly set up." <<std::endl;
        fRecoBranchName="";
        fSimBranchName="";
        fRecoBranchName="";
        fInitialized=false;
     }
}
void createCanvas()
{
     TCanvas *c = (TCanvas*)gROOT->FindObject("c");
     if (c) delete c;
     c = new TCanvas("c","",0,0,1500,700);
     c->Divide(2);
     c->cd(1);
     gPad->SetLeftMargin(0.12);
     c->cd(2);
     gPad->SetLeftMargin(0.12);
     c->cd(1);
}
void setCollimator(const string& tcl, const float& pos)
{
     if (pos==0.) {
        TCut *obj = (TCut*)gROOT->GetList()->FindObject("TCL_F"); gROOT->GetList()->Remove(obj);
        obj = (TCut*)gROOT->GetList()->FindObject("TCL_B"); gROOT->GetList()->Remove(obj);
        TCL_F="";TCL_B=""; TCL_F.SetName("");TCL_B.SetName("");
        setDetectorParameters(detWidth,detHeight,detZPos1,detZPos2,detXPosSigma);
        return;
     }
     TCL_F=Form("%s.ArmF.Xat%s<%f",fSimBranchName.c_str(),tcl.c_str(),pos);
     TCL_B=Form("%s.ArmB.Xat%s<%f",fSimBranchName.c_str(),tcl.c_str(),pos);
     TCL_F.SetName("TCL_F");
     TCL_B.SetName("TCL_B");
     gROOT->GetList()->Add(&TCL_F);
     gROOT->GetList()->Add(&TCL_B);
     setDetectorParameters(detWidth,detHeight,detZPos1,detZPos2,detXPosSigma);
     detF+=TCL_F;
     detB+=TCL_B;
}

void setDetectorParameters(float det_w,float det_h,float det_z1,float det_z2,float det_p, std::string br)
{
// det_p = detector position (in number of sigma of beam)
// det_w = detector width(x)
// det_h = detector height(y)
   if (!fInitialized) setDataFormat();
   if (!fInitialized) {
      std::cout << "Failed to initialize..." << std::endl;
      exit(1);
   }
   if (BeamX.count(1)==0||BeamX.count(2)==0||BeamX.count(3)==0||
       BeamY.count(1)==0||BeamY.count(2)==0||BeamY.count(3)==0||
       BeamX_RMS.count(1)==0||BeamX_RMS.count(2)==0||BeamX_RMS.count(3)==0||
       BeamY_RMS.count(1)==0||BeamY_RMS.count(2)==0||BeamY_RMS.count(3)==0) {
       std::cout << "ERROR: Beam parameters not given..." << std::endl;
       return;
   }
   fSingleArmAcceptance = true;
   //ToFAcceptance = false;
   if (br=="Reco") fTarget=fRecoBranchName;
   else if (br=="Sim") fTarget=fSimBranchName; // make the analysis based on the sim branch

   detXPosSigma=det_p;
   detWidth  = det_w;
   detHeight = det_h;
   detXPos1   = -(BeamX_RMS[1]*detXPosSigma+RPWindowThickness); // detector resides in the negative X
   detXPos2   = -(BeamX_RMS[2]*detXPosSigma+RPWindowThickness); // idem
   detZPos1   = det_z1;
   detZPos2   = det_z2;

   int sig = (detXPos1<0)?-1:1;

   stoppedF=Form("(%s.ArmF.TrkDet1.HasStopped[0]+%s.ArmF.TrkDet2.HasStopped[0])==0",fSimBranchName.c_str(),fSimBranchName.c_str());
   stoppedB=Form("(%s.ArmB.TrkDet1.HasStopped[0]+%s.ArmB.TrkDet2.HasStopped[0])==0",fSimBranchName.c_str(),fSimBranchName.c_str());
   TCut det1x=Form("abs(%s.ArmF.TrkDet1.X[0]-(%f))<%f",fTarget.c_str(),float(sig*(fabs(detXPos1)+det_w/2.)),det_w/2.);
   TCut det1y=Form("abs(%s.ArmF.TrkDet1.Y[0])<%f",fTarget.c_str(),det_h/2);
   TCut det2x=Form("abs(%s.ArmF.TrkDet2.X[0]-(%f))<%f",fTarget.c_str(),float(sig*(fabs(detXPos2)+det_w/2.)),det_w/2.);
   TCut det2y=Form("abs(%s.ArmF.TrkDet2.Y[0])<%f",fTarget.c_str(),det_h/2);
   det1F= det1x&&det1y&&Form("%s.ArmF.TrkDet1.HasStopped[0]==0",fTarget.c_str());
   det2F= det2x&&det2y&&Form("%s.ArmF.TrkDet2.HasStopped[0]==0",fTarget.c_str());
   detF=  det1F&&det2F;
   det1x=Form("abs(%s.ArmB.TrkDet1.X[0]-(%f))<%f",fTarget.c_str(),float(sig*(fabs(detXPos1)+det_w/2)),det_w/2.);
   det1y=Form("abs(%s.ArmB.TrkDet1.Y[0])<%f",fTarget.c_str(),det_h/2);
   det2x=Form("abs(%s.ArmB.TrkDet2.X[0]-(%f))<%f",fTarget.c_str(),float(sig*(fabs(detXPos2)+det_w/2)),det_w/2.);
   det2y=Form("abs(%s.ArmB.TrkDet2.Y[0])<%f",fTarget.c_str(),det_h/2);
   det1B= det1x&&det1y&&Form("%s.ArmB.TrkDet1.HasStopped[0]==0",fTarget.c_str());
   det2B= det2x&&det2y&&Form("%s.ArmB.TrkDet2.HasStopped[0]==0",fTarget.c_str());
   detB=  det1B&&det2B;
   momB= Form("%s.ArmB.momentum[0]>5000&&%s.ArmB.momentum[0]<7000",fTarget.c_str(),fTarget.c_str());
   momF= Form("%s.ArmF.momentum[0]>5000&&%s.ArmF.momentum[0]<7000",fTarget.c_str(),fTarget.c_str());
   detF.SetName("detF");
   detB.SetName("detB");
   det1F.SetName("det1F");
   det2F.SetName("det2F");
   det1B.SetName("det1B");
   det2B.SetName("det2B");
   stoppedF.SetName("stoppedF");
   stoppedB.SetName("stoppedB");
   det1x=Form("abs(%s.ArmF.ToFDet.X[0]-(%f))<%f",fTarget.c_str(),float(sig*(fabs(detXPos1)+det_w/2)),det_w/2.);
   det1y=Form("abs(%s.ArmF.ToFDet.Y[0])<%f",fTarget.c_str(),det_h/2);
   ToFF=det1x&&det1y&&Form("%s.ArmF.ToFDet.HasStopped[0]==0",fTarget.c_str());
   det1x=Form("abs(%s.ArmB.ToFDet.X[0]-(%f))<%f",fTarget.c_str(),float(sig*(fabs(detXPos2)+det_w/2)),det_w/2.);
   det1y=Form("abs(%s.ArmB.ToFDet.Y[0])<%f",fTarget.c_str(),det_h/2);
   ToFB=det1x&&det1y&&Form("%s.ArmB.ToFDet.HasStopped[0]==0",fTarget.c_str());
   momB.SetName("momB");
   momF.SetName("momF");
   TCut *cut = 0;
//
   cut = (TCut*)gROOT->GetList()->FindObject("detF"); if(!cut) gROOT->GetList()->Add(&detF);
   cut = (TCut*)gROOT->GetList()->FindObject("detB"); if(!cut) gROOT->GetList()->Add(&detB);
   cut = (TCut*)gROOT->GetList()->FindObject("momF"); if(!cut) gROOT->GetList()->Add(&momF);
   cut = (TCut*)gROOT->GetList()->FindObject("momB"); if(!cut) gROOT->GetList()->Add(&momB);
   cut = (TCut*)gROOT->GetList()->FindObject("det1F");if(!cut) gROOT->GetList()->Add(&det1F);
   cut = (TCut*)gROOT->GetList()->FindObject("det2F");if(!cut) gROOT->GetList()->Add(&det2F);
   cut = (TCut*)gROOT->GetList()->FindObject("det1B");if(!cut) gROOT->GetList()->Add(&det1B);
   cut = (TCut*)gROOT->GetList()->FindObject("det2B");if(!cut) gROOT->GetList()->Add(&det2B);
   cut = (TCut*)gROOT->GetList()->FindObject("ToFF");if(!cut) gROOT->GetList()->Add(&ToFF);
   cut = (TCut*)gROOT->GetList()->FindObject("ToFB");if(!cut) gROOT->GetList()->Add(&ToFB);
   cut = (TCut*)gROOT->GetList()->FindObject("stoppedF");if(!cut) gROOT->GetList()->Add(&stoppedF);
   cut = (TCut*)gROOT->GetList()->FindObject("stoppedB");if(!cut) gROOT->GetList()->Add(&stoppedB);
}

void drawHitmap(int norm)
{
   TCanvas *c = (TCanvas*)gDirectory->FindObject("c");
   if (c) delete c;
   c = new TCanvas("c","Acceptance",0,0,1000,500);
   c->Divide(2,1);
   c->cd(1);
   drawHitmap("F",norm);
   c->cd(2);
   drawHitmap("B",norm);
   c->cd(0);
}
void drawHitmap(std::string arm,int norm,float xi_min,float xi_max,float xmin,float xmax,float ymin, float ymax)
{
  if (arm!="F"&&arm!="B") {
     std::cout << "arm must be 'F' for forward or 'B' for backward" << std::endl;
     return;
  }

  int nbins = 100;

  gStyle->SetOptStat(0);
  TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
  if (!T||!fInitialized) {
      std::cout << "Environment not properly set up. Quiting." << std::endl;
      return;
  }
  //setZPosition(T);

  float detz = 0;
  if (arm=="F") detz=detZPos1;
  else detz = -detZPos1;

  TVirtualPad *pad = (TVirtualPad*)gROOT->GetSelectedPad();
  if (pad) {pad->SetLogz();}
  else {
     TCanvas *c = new TCanvas("c","Hitmap",0,0,600,600);
     c->SetLogz();
  } 
  
  TCut stopped = *(TCut*)gROOT->GetList()->FindObject(Form("stopped%s",arm.c_str()));
  TCut detgeo = *(TCut*)gROOT->GetList()->FindObject(Form("det1%s",arm.c_str()));
  if (!fSingleArmAcceptance) {
       detgeo = *(TCut*)gROOT->GetList()->FindObject("detF")&&(*(TCut*)gROOT->GetList()->FindObject("detB"));
  }
  TCut tmax   = TCut(Form("%s.Arm%s.t[0]<2",fGenBranchName.c_str(),arm.c_str()));
  TH2F *hxy = new TH2F("hxy","",nbins,xmin,xmax,nbins,ymin,ymax);
  if (xi_min==0&&xi_max==0) {
     T->Draw(Form("%s.Arm%s.TrkDet1.Y[0]:%s.Arm%s.TrkDet1.X[0]>>hxy",fRecoBranchName.c_str(),arm.c_str(),fRecoBranchName.c_str(),arm.c_str()),stopped,"goff");
     hxy->SetTitle(Form("Z_{1}=%3.1f m;Z_{2}=%3.1f m",detZPos1,detZPos2));
  } else {
     T->Draw(Form("%s.Arm%s.TrkDet1.Y[0]:%s.Arm%s.TrkDet1.X[0]>>hxy",fRecoBranchName.c_str(),arm.c_str(),fRecoBranchName.c_str(),arm.c_str()),
             stopped&&Form("%s.Arm%s.xi[0]>%f&&%s.Arm%s.xi[0]<%f",fSimBranchName.c_str(),arm.c_str(),xi_min,fSimBranchName.c_str(),arm.c_str(),xi_max),"goff");
     hxy->SetTitle(Form("Hit map for %3.2f < #xi < %3.2f (z=%+3.1f m)",xi_min,xi_max,detz));
  }
  double generated= 0.;
  double accepted = 0.;
  if (norm) {
     T->Draw(Form("%s.Arm%s.TrkDet1.Y[0]:%s.Arm%s.TrkDet1.X[0]>>hxyTot(%d,%f,%f,%d,%f,%f)",fRecoBranchName.c_str(),arm.c_str(),fRecoBranchName.c_str(),arm.c_str(),
                                                                nbins,xmin,xmax,nbins,ymin,ymax),stopped&&detgeo,"goff");
     TH2F *hxyTot = (TH2F*)gDirectory->Get("hxyTot");
     if (!hxyTot) {
        std::cout<< "Hitmap with accepted protons not found. Quiting." << std::endl;
        return;
     }
     accepted= hxyTot->Integral();
     generated=T->GetEntries();
/*
     for(int i=1;i<=nbins;i++) {
        for(int j=1;j<=nbins;j++) {
           //double gen = hxyTot->GetBinContent(i,j);
           double rec = hxy->GetBinContent(i,j);
           if (gen>0) hxy->SetBinContent(i,j,rec/gen*100.);
        }
     }
*/
  }

  hxy->GetXaxis()->SetTitle("X (mm)");
  hxy->GetYaxis()->SetTitle("Y (mm)");
  hxy->Smooth();
  hxy->Draw("colz");
  //gROOT->GetSelectedPad()->Modified(); gROOT->GetSelectedPad()->Update();
  gPad->Modified(); gPad->Update();
  TPaletteAxis *pa = (TPaletteAxis*)hxy->GetListOfFunctions()->FindObject("palette");
  pa->SetX2NDC(0.93);
  drawDetectorShape();
  if (norm) {
     TText *text = new TText(detXPos1*1.1,-(detHeight/2.)*.9,Form("%2.1f%%",accepted/generated*100));
     text->Draw();
  }
  gROOT->GetSelectedPad()->Modified(); gROOT->GetSelectedPad()->Update();
}

void draw_t_xi_acceptance(int norm, int nbins,double lot,double hit,double loxi,double hixi,std::string fmt)
{
   TCanvas *c = (TCanvas*)gDirectory->FindObject("c");
   if (c) delete c;
   c = new TCanvas("c","Acceptance",0,0,1000,500);
   c->Divide(2,1);
   c->cd(1);
   draw_t_xi_acceptance("F",norm,nbins,lot,hit,loxi,hixi,fmt);
   c->cd(2);
   draw_t_xi_acceptance("B",norm,nbins,lot,hit,loxi,hixi,fmt);
}

void draw_t_xi_acceptance(const std::string& arm, int norm, int nbins,double lot,double hit,double loxi,double hixi,std::string fmt)
{
   gStyle->SetOptStat(kFALSE);
   //gStyle->SetPalette(1);
   if (arm!="F"&&arm!="B") {
      std::cout << "ERROR: wrong arm. It should be 'F' or 'B'" << std::endl;
      return;
   }
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      std::cout << "No tree found. Quiting." << std::endl;
      return;
   }
   //setZPosition(T);
   TVirtualPad *pad = (TVirtualPad*)gROOT->GetSelectedPad();
   if (!pad) {
     TCanvas *c = new TCanvas("c","Hitmap",0,0,600,600);
     c->cd(1);
   } 
   string xi_fmt("");
   string t_fmt("");
   if (fmt.substr(0,3)=="log") {
      xi_fmt="log10";
      if (loxi==0) loxi=0.0001;
      loxi=log10(loxi);hixi=log10(hixi);
   }
   if (fmt=="loglog") {
      if (lot==0) lot=0.0001;
      lot=log10(lot);hit=log10(hit);
      t_fmt="log10";
   }
   
   TH2F* h=(TH2F*)gDirectory->FindObject(Form("hacc%s",arm.c_str())); if (h) delete h;
   //
   string DetType = "";
   TCut XYCut = "";
   TCut XYCutF = "";
   TCut XYCutB = "";
   if (fToFAcceptance) {
      XYCutF = ToFF;
      XYCutB = ToFB;
      DetType = Form("ToF (%2dx%2d@%2.1fmm)(%2.0f#sigma)Z=%3.1fm",int(detWidth),int(detHeight),ToFXPos,detXPosSigma,ToFZPos);
   }
   else {
      XYCutF = detF;
      XYCutB = detB;
      DetType = Form("(%2dx%2d@%2.1fmm)(%2.0f#sigma)",int(detWidth),int(detHeight),detXPos1,detXPosSigma);
   }
   if (arm=="F") XYCut=XYCutF;
   else          XYCut=XYCutB;

   TCut *tcl = (TCut*)gROOT->GetList()->FindObject(Form("TCL_%s",arm.c_str()));
   if (tcl) {
      TString tcl_s(*tcl);
      tcl_s.Remove(0,tcl_s.Index('<')+1);
      double tclpos = tcl_s.Atof();
      TString tcl_id=tcl->GetTitle();
      tcl_id.Remove(0,tcl_id.Index('<')-4);
      tcl_id.Remove(tcl_id.Index('<'));
      DetType+=Form("(%s=%d#sigma)",tcl_id.Data(),int(tclpos));
   }
   TCut stopped = *(TCut*)gROOT->GetList()->FindObject(Form("stopped%s",arm.c_str()));

   T->Draw(Form("%s(%s.Arm%s.t[0]):%s(%s.Arm%s.xi[0])>>hgen%s(%d,%f,%f,%d,%f,%f)",t_fmt.c_str(),fGenBranchName.c_str(),arm.c_str(),xi_fmt.c_str(),fGenBranchName.c_str(),
                                                    arm.c_str(),arm.c_str(),nbins,loxi,hixi,nbins,lot,hit),"","goff");
   if (fSingleArmAcceptance) {
      std::cout << "Acceptance for single arm" << std::endl;
      T->Draw(Form("%s(%s.Arm%s.t[0]):%s(%s.Arm%s.xi[0])>>hrec%s(%d,%f,%f,%d,%f,%f)",t_fmt.c_str(),fGenBranchName.c_str(),arm.c_str(),xi_fmt.c_str(),fGenBranchName.c_str(),
                                                    arm.c_str(),arm.c_str(),nbins,loxi,hixi,nbins,lot,hit),XYCut&&stopped);
   } else {
      std::cout << "Acceptance for double arms" << std::endl;
      T->Draw(Form("%s(%s.Arm%s.t[0]):%s(%s.Arm%s.xi[0])>>hrec%s(%d,%f,%f,%d,%f,%f)",t_fmt.c_str(),fGenBranchName.c_str(),arm.c_str(),xi_fmt.c_str(),fGenBranchName.c_str(),
                                                    arm.c_str(),arm.c_str(),nbins,loxi,hixi,nbins,lot,hit),XYCutF&&stoppedF&&XYCutB&&stoppedB);
   }
   TH2F* hgen = (TH2F*)gDirectory->Get(Form("hgen%s",arm.c_str())); 
   TH2F* hrec = (TH2F*)gDirectory->Get(Form("hrec%s",arm.c_str()));
   if (!hrec||!hgen) {
      std::cout << "Error creating histogram. Quiting." << std::endl;
      return;
   }
   
   TH2F *hacc = new TH2F(Form("hacc%s",arm.c_str()),Form("%s(Z_{1}=%+3.1f m,Z_{2}=%+3.1f m)",DetType.c_str(),detZPos1,detZPos2),nbins,loxi,hixi,nbins,lot,hit);
   
   hacc->GetYaxis()->SetTitleOffset(1.2);hacc->GetXaxis()->SetNdivisions(405);
   if (xi_fmt=="log10") {
      hacc->GetXaxis()->SetTitle("log(#xi)");
      hacc->GetYaxis()->SetTitle("log(t)");
   } else {
      hacc->GetXaxis()->SetTitle("#xi");
      hacc->GetYaxis()->SetTitle("t (GeV/c)^{2}");
   }
   
   for(int i=1;i<=nbins;i++) {
      for(int j=1;j<=nbins;j++) {
         double rec=hrec->GetBinContent(i,j);
         double gen=hgen->GetBinContent(i,j);
         if (gen>0) hacc->SetBinContent(i,j,rec/gen*100.);
      }
   }
   hacc->Smooth();hacc->Smooth();
   hacc->Smooth();hacc->Smooth();
   hacc->SetMinimum(0.);
   hacc->SetMaximum(norm);
   hacc->Draw("colz");
   gPad->Modified(); gPad->Update();

   TPaletteAxis* pa = NULL;
   pa = (TPaletteAxis*)hacc->GetListOfFunctions()->FindObject("palette");
   pa->SetX2NDC(0.93);
   gPad->Modified(); gPad->Update();
}

void draw_totemhitmap(double RPxOffset,double RPyOffset,double xmin,double xmax,double ymin,double ymax)
{
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      std::cout << "Tree not found in current directory. Quiting." << std::endl;
      return;
   }
   setZPosition(T);
   Double_t RPBeamFace = 20.;
   gStyle->SetOptStat(kFALSE);
   gStyle->SetPalette(1);
   TCanvas *c = (TCanvas*)gDirectory->FindObject("c");
   if (c) delete c;
   c = new TCanvas("c","Acceptance",0,0,1000,500);
   c->Divide(2,1);
   c->cd(1);
   TPad* pad = (TPad*)gROOT->FindObject("c_1");
   pad->SetLogz();
   T->Draw(Form("%s.ArmF.ChkPoint.Y:%s.ArmF.ChkPoint.X>>xytotemF(100,%f,%f,100,%f,%f)",fSimBranchName.c_str(),fSimBranchName.c_str(),
                                                                    xmin,xmax,ymin,ymax),det1F&&stoppedF,"colz");
   TH2F* xytotemF = (TH2F*)gDirectory->Get("xytotemF");
   if (!xytotemF) {
      std::cout << "Histogram with hitmap at TOTEM not found. Quiting." <<std::endl;
      return;
   }
   xytotemF->SetTitle(Form("X #times Y at 220m for the accepted at +%3.1fm",detZPos1)); 
   xytotemF->GetXaxis()->SetTitle("X (mm)");
   xytotemF->GetYaxis()->SetTitle("Y (mm)");
   xytotemF->GetYaxis()->SetTitleOffset(1.2);
   xytotemF->Smooth();
   xytotemF->Draw("colz");
   pad->Modified(); pad->Update();
   TPaletteAxis* pa = NULL;
   pa = (TPaletteAxis*)xytotemF->GetListOfFunctions()->FindObject("palette");
   if (pa) pa->SetX2NDC(0.93);
   const int nptv = 4;
   Double_t xplv[nptv] = {xmin,(-RPBeamFace/2<xmin)?xmin:-RPBeamFace/2.,RPBeamFace/2.,(RPBeamFace/2+ymax-RPyOffset>xmax)?xmax:RPBeamFace/2+ymax-RPyOffset};
   Double_t yplv[nptv] = {ymax,RPyOffset,RPyOffset,ymax};
   TPolyLine *ttvu = new TPolyLine(nptv,xplv,yplv);
   ttvu->SetFillStyle(3444);
   ttvu->SetLineColor(kBlack);
   ttvu->SetFillColor(38);
   ttvu->SetLineWidth(2);
   for(int i=0;i<nptv;i++) yplv[i]*=-1;
   TPolyLine *ttvd = new TPolyLine(nptv,xplv,yplv);
   ttvd->SetFillStyle(3444);
   ttvd->SetLineColor(kBlack);
   ttvd->SetFillColor(38);
   ttvd->SetLineWidth(2);
   Double_t xplh[6] = {xmax,RPxOffset+ymax-RPBeamFace/2,RPxOffset,RPxOffset,RPxOffset+ymax-RPBeamFace/2.,xmax};
   Double_t yplh[6] = {ymax,ymax,RPBeamFace/2.,-RPBeamFace/2.,ymin,ymin};
   TPolyLine *tth = new TPolyLine(6,xplh,yplh);
   tth->SetFillStyle(3444);
   tth->SetLineColor(kBlack);
   tth->SetFillColor(38);
   tth->SetLineWidth(2);
   ttvu->Draw("f"); ttvu->Draw();
   ttvd->Draw("f"); ttvd->Draw();
   tth->Draw("f"); tth->Draw();
   pad->Modified(); pad->Update();
   
   c->cd(2);
   pad = (TPad*)gROOT->FindObject("c_2");
   pad->SetLogz();
   T->Draw(Form("%s.ArmB.ChkPoint.Y:%s.ArmB.ChkPoint.X>>xytotemB(100,%f,%f,100,%f,%f)",fSimBranchName.c_str(),fSimBranchName.c_str(),
                                                                    xmin,xmax,ymin,ymax),det1B&&stoppedB,"colz");
   TH2F* xytotemB = (TH2F*)gDirectory->Get("xytotemB");
   if (!xytotemB) {
      std::cout << "Histogram with TOTEM hitmap not found. Quiting."<<std::endl;
      return;
   }
   xytotemB->SetMaximum(1000);
   xytotemB->SetTitle(Form("X #times Y at -220m for the accepted at -%3.1fm",detZPos1)); 
   xytotemB->GetXaxis()->SetTitle("X (mm)");
   xytotemB->GetYaxis()->SetTitleOffset(1.2);
   xytotemB->GetYaxis()->SetTitle("Y (mm)");
   xytotemB->Smooth();
   pad->Modified(); pad->Update();
   pa = (TPaletteAxis*)xytotemB->GetListOfFunctions()->FindObject("palette");
   if (pa) pa->SetX2NDC(0.93);
   ttvu->Draw("f"); ttvu->Draw();
   ttvd->Draw("f"); ttvd->Draw();
   tth->Draw("f"); tth->Draw();
   pad->Modified(); pad->Update();
   c->Update();
}
void draw_t_vs_mass_acceptance(const char* arm, double norm,int nbins, double mmin, double mmax,double tmin, double tmax,double ecms)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(1);
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      std::cout << "No tree found. Quiting." << std::endl;
      return;
   }
   setZPosition(T);
   T->Draw(Form("%s.Arm%s.t:sqrt(%s.ArmB.xi[0]*%s.ArmF.xi[0])*%f>>mxt_gen%s(%d,%f,%f,%d,%f,%f)",
                fGenBranchName.c_str(),arm,fGenBranchName.c_str(),arm,
                ecms,fGenBranchName.c_str(),nbins,mmin,mmax,nbins,tmin,tmax));
   T->Draw(Form("%s.Arm%s.t:sqrt(%s.ArmB.xi[0]*%s.ArmF.xi[0])*%f>>mxt_reco%s(%d,%f,%f,%d,%f,%f)",
                 fGenBranchName.c_str(),arm,fGenBranchName.c_str(),arm,
                 ecms,fGenBranchName.c_str(),nbins,mmin,mmax,nbins,tmin,tmax),detF&&detB&&stoppedF&&stoppedB);
   TH2F *mxt_gen = (TH2F*)gDirectory->Get(Form("mxt_gen%s",arm));
   if (!mxt_gen) {
      std::cout<< "Central mass with generated protons not found. Quiting." << std::endl;
      return;
   }
   TH2F *mxt_reco = (TH2F*)gDirectory->Get(Form("mxt_reco%s",arm));
   if (!mxt_reco) {
      std::cout<< "Central mass with reconstructed protons not found. Quiting." << std::endl;
      return;
   }
   TH2F *mxt_acc = (TH2F*)gDirectory->Get(Form("mxt_acc%s",arm));
   if (mxt_acc) delete mxt_acc;
   mxt_acc = new TH2F(Form("mxt_acc%s",arm),"",nbins,mmin,mmax,nbins,tmin,tmax);

   for(int i=1;i<=nbins;i++) {
      for(int j=1;j<=nbins;j++) {
         double mgen = mxt_gen->GetBinContent(i,j);
         double mreco =mxt_reco->GetBinContent(i,j);
         if (mgen>0) mxt_acc->SetBinContent(i,j,mreco/mgen*100);
      }
   }
   mxt_acc->Smooth();
   mxt_acc->Draw("colz");
   if (norm<mxt_acc->GetMaximum()) std::cout << "WARNING: setting the histogram to a maximum smaller than its maximum" << std::endl;
   mxt_acc->SetMaximum(norm);
   mxt_acc->Draw("colz");
   if (*arm=='F')
      mxt_acc->SetTitle(Form("t #times M_{X} acceptance at Z_{1}=%3.1fm and Z_{2}=%3.1fm;M_{X} (GeV/c^{2});t (GeV/c)^{2}",detZPos1,detZPos2));
   else
      mxt_acc->SetTitle(Form("t #times M_{X} acceptance at Z_{1}=%3.1fm and Z_{2}=%3.1fm;M_{X} (GeV/c^{2});t (GeV/c)^{2}",-detZPos1,-detZPos2));

   gPad->Modified();gPad->Update();
   TPaletteAxis* pa = NULL;
   pa = (TPaletteAxis*)mxt_acc->GetListOfFunctions()->FindObject("palette");
   pa->SetX2NDC(0.93);
   gPad->Modified();gPad->Update();
}
void do_certification(double min_t,double max_t, double min_xi, double max_xi)
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.2);
   gStyle->SetPalette(1);
   
   TCanvas *t = (TCanvas*)gDirectory->FindObject("t");
   if (t) delete t;
   t = new TCanvas("t","",0,0,1000,500);
   t->Divide(2,1);
   t->cd(1);
   TPad* pad = (TPad*)gROOT->FindObject("t_1");
   pad->SetLogz();
   draw_t_correlation("F",min_t,max_t);
   t->cd(2);
   pad = (TPad*)gROOT->FindObject("t_2");
   pad->SetLogz();
   draw_t_correlation("B",min_t,max_t);
//
// xi correlations
   TCanvas *x= (TCanvas*)gDirectory->FindObject("x");
   if (x) delete x;
   x = new TCanvas("x","",1100,0,1000,500);
   x->Divide(2,1);
   x->cd(1);
   pad = (TPad*)gROOT->FindObject("x_1");
   pad->SetLogz();
   draw_xi_correlation("F",min_xi,max_xi);
   x->cd(2);
   pad = (TPad*)gROOT->FindObject("x_2");
   pad->SetLogz();
   draw_xi_correlation("B",min_xi,max_xi);
}

void draw_t_correlation(const char* arm,double min_t,double max_t)
{
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      std::cout << "No tree found. Quiting." << std::endl;
      return;
   }
   setZPosition(T);
   double det_z = 0.;
   TCut det = "";
   TCut mom = "";
   if (*arm=='F') {
      det_z = detZPos1;
      det = detF;
      mom = momF;
   }
   if (*arm=='B') {
      det_z = -detZPos1;
      det = detB;
      mom = momB;
   }
   T->Draw(Form("%s.Arm%s.t:%s.Arm%s.t>>ht%s(100,%f,%f,100,%f,%f)",fTarget.c_str(),arm,fGenBranchName.c_str(),
                                          arm,arm,min_t,max_t,min_t,max_t),det&&mom,"colz");
   TH2F *ht = (TH2F*)gDirectory->FindObject(Form("ht%s",arm));
   ht->Smooth();
   ht->SetTitle(Form("Z = %+3.1f m;t_{gen} (GeV^{2});t_{reco} (GeV^{2})",det_z));
   ht->Draw("colz");
   gPad->Modified();gPad->Update();
   TPaletteAxis* pa = NULL;
   pa = (TPaletteAxis*)ht->GetListOfFunctions()->FindObject("palette");
   pa->SetX2NDC(0.93);
   gPad->Modified();gPad->Update();
}
void draw_xi_correlation(const char* arm,double min_xi,double max_xi)
{
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      std::cout << "No tree found. Quiting." << std::endl;
      return;
   }
   setZPosition(T);
   double det_z = 0.;
   TCut det = "";
   TCut mom = "";
   if (*arm=='F') {
      det_z = detZPos1;
      det = detF;
      mom = momF;
   }
   if (*arm=='B') {
      det_z = -detZPos1;
      det = detB;
      mom = momB;
   }
   T->Draw(Form("%s.Arm%s.xi:%s.Arm%s.xi>>hxi%s(100,%f,%f,100,%f,%f)",fTarget.c_str(),arm,fGenBranchName.c_str(),
                                  arm,arm,min_xi,max_xi,min_xi,max_xi),det&&mom,"colz");
   TH2F *hxi = (TH2F*)gDirectory->FindObject(Form("hxi%s",arm));
   hxi->Smooth();

   hxi->SetTitle(Form("Z = %+3.1f m;#xi_{gen};#xi_{reco}",det_z));
   hxi->Draw("colz");
   gPad->Modified();gPad->Update();
   TPaletteAxis* pa = NULL;
   pa = (TPaletteAxis*)hxi->GetListOfFunctions()->FindObject("palette");
   pa->SetX2NDC(0.93);
   gPad->Modified();gPad->Update();
}
void draw_y_vs_mass_acceptance(double norm,int nbins, double mmin, double mmax, double ymin, double ymax,double ecms)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(1);
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      std::cout << "No tree found. Quiting." << std::endl;
      return;
   }
   setZPosition(T);
   T->Draw(Form("abs(log(%s.ArmF.xi[0]/%s.ArmB.xi[0])/2.0):sqrt(%s.ArmF.xi[0]*%s.ArmB.xi[0])*%f>>ymxgen(%d,%f,%f,%d,%f,%f)",
                        fGenBranchName.c_str(),fGenBranchName.c_str(),fGenBranchName.c_str(),fGenBranchName.c_str(),
                                                                     ecms,nbins,mmin,mmax,nbins,ymin,ymax),"","colz");
   T->Draw(Form("abs(log(%s.ArmF.xi[0]/%s.ArmB.xi[0])/2.0):sqrt(%s.ArmF.xi[0]*%s.ArmB.xi[0])*%f>>ymxrec(%d,%f,%f,%d,%f,%f)",
                        fGenBranchName.c_str(),fGenBranchName.c_str(),fGenBranchName.c_str(),fGenBranchName.c_str(),
                                                                     ecms,nbins,mmin,mmax,nbins,ymin,ymax),
                                                                           detF&&detB&&momF&&momB,"colz");
   TH2F* ymxgen = (TH2F*)gDirectory->FindObject("ymxgen");
   TH2F* ymxrec = (TH2F*)gDirectory->FindObject("ymxrec");
   TH2F* ymxacc = (TH2F*)gDirectory->FindObject("ymxacc");
   if (ymxacc) delete ymxacc;
   ymxacc = new TH2F("ymxacc","",nbins,mmin,mmax,nbins,ymin,ymax);
   for(int i=1;i<=nbins;i++) {
      for(int j=1;j<=nbins;j++) {
         double gen=ymxgen->GetBinContent(i,j);
         double rec=ymxrec->GetBinContent(i,j);
         if (gen>0) ymxacc->SetBinContent(i,j,double(rec/gen)*100.);
      }
   }
   ymxacc->Smooth();
   ymxacc->Smooth();
   ymxacc->SetStats(0);
   if (norm<ymxacc->GetMaximum()) std::cout << "WARNING: setting the histogram to a maximum smaller than its maximum" << std::endl;
   ymxacc->SetMaximum(norm);
   ymxacc->SetTitle(Form("rapidity #times M_{X} acceptance at Z_{1}=%3.1f m and Z_{2}=%3.1f m;M_{X} (GeV/c^{2});y",detZPos1,detZPos2));
   ymxacc->Draw("colz");
   gPad->Modified();gPad->Update();
   TPaletteAxis *pa = (TPaletteAxis*)ymxacc->GetListOfFunctions()->FindObject("palette");
   pa->SetX2NDC(0.93);
   gPad->Modified();gPad->Update();
}
void draw_mass_acceptance(TCut cut,const char* opt, double ecms, double norm, int nbins, double mmin, double mmax)
{
   gStyle->SetOptStat(0);
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      std::cout << "No tree found. Quiting." << std::endl;
      return;
   }
   setZPosition(T);
   T->Draw(Form("sqrt(%s.ArmF.xi[0]*%s.ArmB.xi[0])*%f>>hgen(%d,%f,%f)",fGenBranchName.c_str(),fGenBranchName.c_str(),ecms,
                                                    nbins,mmin,mmax),cut,"goff");
   T->Draw(Form("sqrt(%s.ArmF.xi[0]*%s.ArmB.xi[0])*%f>>hrec(%d,%f,%f)",fGenBranchName.c_str(),fGenBranchName.c_str(),ecms,
                                                    nbins,mmin,mmax),detF&&detB&&cut,"goff");
   TH1F* hgen = (TH1F*)gDirectory->FindObject("hgen");
   TH1F* hrec = (TH1F*)gDirectory->FindObject("hrec");
   TH1F* hacc = (TH1F*)gDirectory->FindObject("hacc");
   if (hacc) delete hacc;
   hacc = new TH1F(*hrec);
   hacc->SetName("hacc");
   hacc->Divide(hgen);
   hacc->Scale(100);
   hacc->SetMaximum(norm);
   hacc->SetStats(0);
   hacc->SetTitle(Form("Mass acceptance (%2dx%2dmm @ %2.1f#sigma), Z_{1}=#pm%3.1fm and Z_{2}=#pm%3.1fm;M_{X} (GeV/c^{2});acceptance(%%)",
                                      int(detWidth),int(detHeight),detXPosSigma,detZPos1,detZPos2));
   if (*opt=='H'||*opt=='h') hacc->Draw();
   else {
      hacc->Smooth();
      hacc->Draw("C");
   }
}
void draw_centralmass(double ecms, int nbins, double mmin, double mmax)
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit();
   //gStyle->SetFitFormat("4.3g");
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      std::cout << "No tree found. Quiting." << std::endl;
      return;
   }
   setZPosition(T);
   if (fRestrictedGen) {
      T->Draw(Form("sqrt(%s.ArmF.xi[0]*%s.ArmB.xi[0])*%f>>mx_gen(%d,%f,%f)",fGenBranchName.c_str(),fGenBranchName.c_str(),ecms,nbins,mmin,mmax),
                                                                   Form("(%s.ArmF.TrkDet1.HasStopped[0]+%s.ArmB.TrkDet1.HasStopped[0])==0",
                                                                   fSimBranchName.c_str(),fSimBranchName.c_str()),"goff");
   } else {
     T->Draw(Form("sqrt(%s.ArmF.xi[0]*%s.ArmB.xi[0])*%f>>mx_gen(%d,%f,%f)",fGenBranchName.c_str(),fGenBranchName.c_str(),ecms,nbins,mmin,mmax));
   }
   T->Draw(Form("sqrt(%s.ArmF.xi[0]*%s.ArmB.xi[0])*%f>>mx_sim(%d,%f,%f)",fSimBranchName.c_str(),fSimBranchName.c_str(),ecms,nbins,mmin,mmax),
                                                                   Form("(%s.ArmF.TrkDet1.HasStopped[0]+%s.ArmB.TrkDet1.HasStopped[0])==0",
                                                                   fSimBranchName.c_str(),fSimBranchName.c_str()),"goff");
   T->Draw(Form("sqrt(%s.ArmF.xi[0]*%s.ArmB.xi[0])*%f>>mx_reco(%d,%f,%f)",fRecoBranchName.c_str(),fRecoBranchName.c_str(),ecms,
                                                                   nbins,mmin,mmax),detF&&momF&&detB&&momB,"goff");
   TH1F *mx_gen = (TH1F*)gDirectory->FindObject("mx_gen");
   TH1F *mx_sim = (TH1F*)gDirectory->FindObject("mx_sim");
   TH1F *mx_reco= (TH1F*)gDirectory->FindObject("mx_reco");
   mx_gen->SetTitle("M_{X} = #sqrt{#xi_{1}#xi_{2}}#sqrt{s};M_{X} (GeV/c^{2});");
   mx_sim->SetTitle("M_{X} = #sqrt{#xi_{1}#xi_{2}}#sqrt{s};M_{X} (GeV/c^{2});");
   mx_reco->SetTitle("M_{X} = #sqrt{#xi_{1}#xi_{2}}#sqrt{s};M_{X} (GeV/c^{2});");
   mx_reco->SetLineWidth(2); mx_gen->SetLineWidth(2); mx_sim->SetLineWidth(2);
   mx_reco->SetLineColor(kRed); mx_gen->SetLineColor(kBlack); mx_sim->SetLineColor(kBlue);
   if (mx_gen->GetEffectiveEntries()>0) mx_gen->Scale(1./mx_gen->GetMaximum());
   if (mx_sim->GetEffectiveEntries()>0) mx_sim->Scale(1./mx_sim->GetMaximum());
   if (mx_reco->GetEffectiveEntries()>0)mx_reco->Scale(1./mx_reco->GetMaximum());
   mx_gen->Draw();mx_sim->Draw("same");mx_reco->Draw("same");
   TLegend *leg = new TLegend(0.50,0.75,0.9,0.9);
   leg->SetFillStyle(0);
   if (fRestrictedGen) leg->AddEntry(mx_gen,"M_{X} Generated at IP (restricted)","l");
   else leg->AddEntry(mx_gen,"M_{X} Generated at CMS","l");
   leg->AddEntry(mx_sim,Form("M_{X} propagated to PPS (Z=%3.1f m)",detZPos1),"l");
   leg->AddEntry(mx_reco,Form("M_{X} reco. at PPS (%2.0f#times%2.0f@%2.1f mm)(%2.0f#sigma)",detWidth,detHeight,detXPos1,detXPosSigma),"l");
   leg->Draw();
   //mx_gen->Fit("gaus");
   //mx_sim->Fit("gaus");
   //mx_reco->Fit("gaus","","sames");
   //mx_reco->GetFunction("gaus")->SetLineColor(kRed);
   //gPad->Modified(); gPad->Update();
   //TPaveStats* stats = (TPaveStats*)mx_sim->FindObject("stats");
   //float stat_offset = 0.07;
   //stats->SetY2NDC(stats->GetY2NDC()-stat_offset);
   //stats->SetY1NDC(stats->GetY1NDC()-stat_offset);
   //stats = (TPaveStats*)mx_reco->FindObject("stats");
   //float boxy1=stats->GetY1NDC();
   //float boxy2=stats->GetY2NDC();
   //stats->SetY1NDC(boxy1-(boxy2-boxy1)-stat_offset);stats->SetY2NDC(boxy1-stat_offset);
   //stats->SetTextColor(kBlue);
   //mx_gen->Draw();mx_sim->Draw("same");mx_reco->Draw("same");
}
void draw_mass_resolution(double mass,double ecms, const char* br,const char* cut, int nbins, double d_mass)
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit();
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      std::cout << "No tree found. Quiting." << std::endl;
      return;
   }
   setZPosition(T);
   TString m_1 = Form("sqrt(%s.ArmF.xi[0]*%s.ArmB.xi[0])*%f",fTarget.c_str(),fTarget.c_str(),ecms);
   TString m_g = Form("sqrt(%s.ArmF.xi[0]*%s.ArmB.xi[0])*%f",fGenBranchName.c_str(),fGenBranchName.c_str(),ecms);

   if (mass>0) {
      T->Draw(Form("(%s-%s)>>mx_res(%d,%f,%f)",m_1.Data(),m_g.Data(),nbins,-d_mass,d_mass),detF&&momF&&detB&&momB&&cut,"goff");
   } else {
      T->Draw(Form("(%s-%s)/(%s)>>mx_res(%d,%f,%f)",m_1.Data(),m_g.Data(),m_g.Data(),nbins,-d_mass,d_mass),
              detF&&momF&&detB&&momB&&cut,"goff");
   } 
   TH1F *mx_res = (TH1F*)gDirectory->FindObject("mx_res");
   mx_res->SetLineWidth(2);mx_res->SetLineColor(kRed);
   mx_res->SetTitle(Form("M_{X} resolution (%2.0f#times%2.0f@%2.1f mm)(%2.0f#sigma) at Z_{1}=%3.1f m,Z_{2}=%3.1f m",
                    detWidth,detHeight,detXPos1,detXPosSigma,detZPos1,detZPos2));
   mx_res->GetYaxis()->SetTitle("a.u.");

   if (mass>0) mx_res->GetXaxis()->SetTitle(Form("M_{X}^{%s}-M_{X}^{Gen} (GeV/c^{2})",br));
   else        mx_res->GetXaxis()->SetTitle(Form("(M_{X}^{%s}-M_{X}^{Gen})/M_{X}^{Gen}",br));

   mx_res->Scale(1./mx_res->GetMaximum());
   mx_res->Fit("gaus");
   mx_res->GetFunction("gaus")->SetLineColor(kBlue);
   gPad->Modified(); gPad->Update();
   TPaveStats* stats = (TPaveStats*)mx_res->FindObject("stats");
   float stat_offset = 0.07;
   stats->SetY2NDC(stats->GetY2NDC()-stat_offset);
   stats->SetY1NDC(stats->GetY1NDC()-stat_offset);
   if (mass>0) {
      TLatex* tmx = new TLatex(-d_mass*0.9,0.95,Form("M_{X}=%4.0f GeV/c^{2}",mass));
      tmx->Draw();
   }
   gPad->Modified();gPad->Update();
}
void setZPosition(TTree* t)
{
     //if (fTreeName!="T") return;
     if (detZPos1>0&&detZPos2>0) return;
     detZPos1 = 0.;
     detZPos2 = 0.;
     if (!t) return;
     std::vector<double> det1Z;
     std::vector<double> det2Z;
     t->SetMakeClass(1);
     t->SetBranchAddress(Form("%s.ArmF.TrkDet1.Z",fSimBranchName.c_str()),&det1Z);
     t->SetBranchAddress(Form("%s.ArmF.TrkDet2.Z",fSimBranchName.c_str()),&det2Z);
     t->GetEntry(0);
     if (!det1Z.size()||!det1Z.size()) return;
     detZPos1 = det1Z.at(0);
     detZPos2 = det2Z.at(0);
     t->ResetBranchAddresses();
     t->SetMakeClass(0);
}
void draw_t_xi_rejection(int norm , int nbins,double lot,double hit,double loxi,double hixi,string fmt)
{
   gStyle->SetOptStat(kFALSE);
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      std::cout << "No tree found. Quiting." << std::endl;
      return;
   }
   setZPosition(T);
   TCanvas *c = (TCanvas*)gDirectory->FindObject("c");
   if (c) delete c;
   c = new TCanvas("c","Acceptance",0,0,1000,500);
   c->Divide(2,1);
   TH2F* h=(TH2F*)gDirectory->FindObject("haccF"); if (h) delete h;
   h=(TH2F*)gDirectory->FindObject("haccB"); if (h) delete h;
//
   T->Draw(Form("%s(%s.ArmF.t[0]):%s(%s.ArmF.xi[0])>>hgenF(%d,%f,%f,%d,%f,%f)",fmt.c_str(),fGenBranchName.c_str(),fmt.c_str(),fGenBranchName.c_str(),
                                                   nbins,loxi,hixi,nbins,lot,hit),PPSanalyzer::det1F,"goff");
   T->Draw(Form("%s(%s.ArmB.t[0]):%s(%s.ArmB.xi[0])>>hgenB(%d,%f,%f,%d,%f,%f)",fmt.c_str(),fGenBranchName.c_str(),fmt.c_str(),fGenBranchName.c_str(),
                                                   nbins,loxi,hixi,nbins,lot,hit),PPSanalyzer::det1B,"goff");
//
   T->Draw(Form("%s(%s.ArmF.t[0]):%s(%s.ArmF.xi[0])>>hrecF(%d,%f,%f,%d,%f,%f)",fmt.c_str(),fGenBranchName.c_str(),fmt.c_str(),fGenBranchName.c_str(),
                                                   nbins,loxi,hixi,nbins,lot,hit),PPSanalyzer::det1F&&(!PPSanalyzer::det2F),"goff");
   T->Draw(Form("%s(%s.ArmB.t[0]):%s(%s.ArmB.xi[0])>>hrecB(%d,%f,%f,%d,%f,%f)",fmt.c_str(),fGenBranchName.c_str(),fmt.c_str(),fGenBranchName.c_str(),
                                                   nbins,loxi,hixi,nbins,lot,hit),PPSanalyzer::det1B&&(!PPSanalyzer::det2B),"goff");

   TH2F* hgenF = (TH2F*)gDirectory->Get("hgenF"); 
   TH2F* hrecF = (TH2F*)gDirectory->Get("hrecF");
   TH2F* hgenB = (TH2F*)gDirectory->Get("hgenB");
   TH2F* hrecB = (TH2F*)gDirectory->Get("hrecB");
   if (!hrecF||!hgenF||!hrecB||!hgenB) {
      std::cout << "Error creating histogram. Quiting." << std::endl;
      return;
   }
   
   TH2F *haccF = new TH2F("haccF",Form("Det_{2} Rejection rate (%2.0f#times%2.0f@%2.1f mm)(%2.0f#sigma) at Z_{1}=%3.1fm and Z_{2}=%3.1fm;#xi;|t| (GeV/c)^{2}",
                     detWidth,detHeight,detXPos1,detXPosSigma,detZPos1,detZPos2),nbins,loxi,hixi,nbins,lot,hit);
   TH2F *haccB = new TH2F("haccB",Form("Det_{2} Rejection rate (%2.0f#times%2.0f@%2.1f mm)(%2.0f#sigma) at Z_{1}=%3.1fm and Z_{2}=%3.1fm;#xi;|t| (GeV/c)^{2}",
                     detWidth,detHeight,(detXPos1),detXPosSigma,-detZPos1,-detZPos2),nbins,loxi,hixi,nbins,lot,hit);

   haccF->GetYaxis()->SetTitleOffset(1.2);haccF->GetXaxis()->SetNdivisions(405);
   haccB->GetYaxis()->SetTitleOffset(1.2);haccB->GetXaxis()->SetNdivisions(405);
   for(int i=1;i<=nbins;i++) {
      for(int j=1;j<=nbins;j++) {
         double gen = hgenF->GetBinContent(i,j);
         double rec = hrecF->GetBinContent(i,j);
         if (gen>0) haccF->SetBinContent(i,j,rec/gen*100);
         else haccF->SetBinContent(i,j,0.);
         gen = hgenB->GetBinContent(i,j);
         rec = hrecB->GetBinContent(i,j);
         if (gen>0) haccB->SetBinContent(i,j,rec/gen*100);
         else haccB->SetBinContent(i,j,0.);
      }
   }
   haccF->Smooth();
   haccB->Smooth();
   haccF->SetMinimum(0.); haccB->SetMinimum(0.);
   haccF->SetMaximum(norm); haccB->SetMaximum(norm);
   c->cd(1); haccF->Draw("colz");
   gPad->Modified(); gPad->Update();

   TPaletteAxis* pa = NULL;
   pa = (TPaletteAxis*)haccF->GetListOfFunctions()->FindObject("palette");
   pa->SetX2NDC(0.93);
   gPad->Modified(); gPad->Update();
   //
   c->cd(2); haccB->Draw("colz");
   gPad->Modified(); gPad->Update();
   pa = (TPaletteAxis*)haccB->GetListOfFunctions()->FindObject("palette");
   pa->SetX2NDC(0.93);
   gPad->Modified(); gPad->Update();
}
void drawDetectorShape()
{
     int sig = (detXPos1<0)?-1:1;
     float xmin = min(detXPos1,detXPos1+sig*detWidth);
     float xmax = max(detXPos1,detXPos1+sig*detWidth);
     TBox *pps = new TBox(xmin,-detHeight/2,xmax,detHeight/2.);
     pps->SetLineWidth(3);
     pps->SetFillStyle(0);
     pps->SetLineColor(kBlack);
     float xtext = detXPos1+sig*detWidth;
     TText *tx = NULL;
     if (sig<0) tx = new TText(xtext+1,-detHeight/2-1.7,Form("%2.1fmm",detWidth));
     else       tx = new TText(xtext-detWidth/3,-detHeight/2-1.7,Form("%2.1fmm",detWidth));
     TText *ty = new TText(xtext+sig,-detHeight/2+0.5,Form("%2.1fmm",detHeight));
     ty->SetTextAngle(90);
     tx->Draw();ty->Draw();
     pps->Draw();
}
void drawBeamEnvelope(double sig)
{
     TEllipse *beamE = new TEllipse(BeamX[1],BeamY[1],sig*BeamX_RMS[1],sig*BeamY_RMS[1],90,270);
     beamE->SetLineWidth(2);
     beamE->SetFillStyle(0);
     beamE->SetNoEdges();
     beamE->Draw();
}
void setToFGeometry(std::string geotype)
{
     ToFCellRow.clear();
     ToFCellColumn.clear();
     if (geotype=="Quartic") {
        ToFGeometryType = "Quartic";
        setQuartic();
     } else if (geotype=="Diamond") {
        ToFGeometryType = "Diamond";
        setDiamond();
     } else {
        std::cout << "Unknown ToF geometry." << std::endl;
     }
}
void setQuartic()
{
// the vertical positions starts from the negative(bottom) to the positive(top) corner
     for(int i=0;i<NYCell;i++) {
// vector index points to the row number from below
        double y1=QuarticCellHeight*(i-NYCell/2.);
        double y2=QuarticCellHeight*(i-NYCell/2.+1);
        ToFCellRow.push_back(pair<double,double>(y1,y2));
     }
// vector index points to the column number
     for(int i=0;i<NXCell;i++) {
        double x1 = ToFXPos-(QuarticCellWidth*i);
        double x2 = ToFXPos-(QuarticCellWidth*(i+1));
        ToFCellColumn.push_back(pair<double,double>(x1,x2));
     }
}
void setDiamond()
{
// vector index points do the row number (3 - bottom/center/top)
    ToFCellRow.push_back(pair<double,double>(-3*DiamondCellHeight/2.,-DiamondCellHeight/2.));
    ToFCellRow.push_back(pair<double,double>(-DiamondCellHeight/2.,DiamondCellHeight/2));
    ToFCellRow.push_back(pair<double,double>(DiamondCellHeight/2.,3*DiamondCellHeight/2));
// define the cell widths (the columns)
    ToFCellColumn.push_back(pair<double,double>(ToFXPos,ToFXPos-DiamondLowerCellW));
// vector index points to the column (variable width) in the central row
    double dx=0;
    for(unsigned int i=0;i<DiamondCentralCellW.size();i++) {
       double x1 = ToFXPos-dx;
       double x2 = ToFXPos-(dx+DiamondCentralCellW.at(i));
       dx+=DiamondCentralCellW.at(i);
       ToFCellColumn.push_back(pair<double,double>(x1,x2));
    }
// vector with with the top and bottom row (1 bigger cell each)
    ToFCellColumn.push_back(pair<double,double>(ToFXPos,ToFXPos-DiamondUpperCellW));
}
int ToFCellId(double x, double y)
{
   int y_idx,x_idx;
// first, get the row number
   unsigned int i;
   unsigned int start_idx=0;
   unsigned int end_idx=ToFCellRow.size();
   for(i=0;i<ToFCellRow.size();i++){
      if (y>=ToFCellRow.at(i).first&&y<=ToFCellRow.at(i).second) break;
   } 
   if (i>=ToFCellRow.size()) return 0;
   y_idx = i+1;
   if (ToFGeometryType=="Quartic") {
      start_idx=0;end_idx=ToFCellColumn.size();
   }
   if (ToFGeometryType=="Diamond") {
      switch (y_idx) {
             case 1: start_idx=0; end_idx=1;break;
             case 3: start_idx=DiamondCentralCellW.size()+1;end_idx=ToFCellColumn.size();break;
             case 2: start_idx=1;end_idx=1+DiamondCentralCellW.size();break;
             default: cout  << "ERROR: unknown ToF row index" << endl;return 0;
      }
   }
   for(i=start_idx;i<end_idx;i++) {
      if (x<=ToFCellColumn.at(i).first&&x>ToFCellColumn.at(i).second) break;
   }
   if (i>=end_idx) return 0;
   x_idx=i+1-start_idx;
   return 100*y_idx+x_idx;
}
double ToFXFromCell(int celid)
{
   if(!isValidCellId(celid)) return 0;
   unsigned int y_idx=int(celid/100);
   unsigned int x_idx=celid-y_idx*100;
   double xhit = 0;
   if (ToFGeometryType=="Quartic") {
      xhit = (ToFCellColumn.at(x_idx-1).first+ToFCellColumn.at(x_idx-1).second)/2.0;
   } else if (ToFGeometryType=="Diamond") {
     switch (y_idx) {
            case 1: xhit = (ToFCellColumn.at(x_idx-1).first+ToFCellColumn.at(x_idx-1).second)/2.0;
                    break;
            case 3: xhit = (ToFCellColumn.at(DiamondCentralCellW.size()+1).first+ToFCellColumn.at(DiamondCentralCellW.size()+1).second)/2.;
                    break;
            case 2: xhit = (ToFCellColumn.at(x_idx).first+ToFCellColumn.at(x_idx).second)/2.;
                    break;
            default: cout  << "ERROR: unknown ToF row index" << endl;return 0;
     }
   } else {
     cout << "ERROR: Unknown ToF geometry." << endl;
     return 0;
   }
   return xhit;
}
double ToFYFromCell(int celid)
{
   if (!isValidCellId(celid)) return 0;
   unsigned int y_idx=int(celid/100);
   double yhit = (ToFCellRow.at(y_idx-1).first+ToFCellRow.at(y_idx-1).second)/2.0;
   return yhit;
}
bool isValidCellId(int celid)
{
     int y_idx=int(celid/100);
     int x_idx=celid-y_idx*100;
     if (ToFGeometryType=="Quartic") {
        if (y_idx<1||y_idx>NYCell) return false;
        if (x_idx<1||x_idx>NXCell) return false;
        return true;
     }
     if (ToFGeometryType=="Diamond") {
        if (y_idx<1||y_idx>DiamondNYCell) return false;
        if (x_idx<1) return false;
        if (y_idx==1&&x_idx>NXCellLowerRow) return false;
        if (y_idx==3&&x_idx>NXCellUpperRow) return false;
        if (y_idx==2&&x_idx>NXCellCentralRow) return false;
        return true;
     }
     cout << "ERROR: Unknown ToF geometry."<< endl;
     return false;
}
void drawToFOccupancy(const string& arm, double xmin, double xmax, double ymin, double ymax)
{
     if (ToFGeometryType=="Quartic") drawQuarticOccupancy(arm,xmin,xmax,ymin,ymax);
     if (ToFGeometryType=="Diamond") drawDiamondOccupancy(arm,xmin,xmax,ymin,ymax);
}
void drawDiamondOccupancy(const string& arm, double xmin, double xmax, double ymin, double ymax)
{
   double xbin[NXCellCentralRow+1+2]; // make room for 2 adicional bins in order to show it in the optimal scale
   int i_skip=0;
   int nbin=NXCellCentralRow;
   int first_tofbin=0;
   if (xmax>ToFXPos) {
      xbin[NXCellCentralRow+2] = xmax;
      xbin[NXCellCentralRow+1] = ToFXPos;
      first_tofbin=1;
      i_skip++;
      nbin++;
   }
   else {
      xbin[NXCellCentralRow]=ToFXPos;
   }
   int i=0;
   for(i=NXCellCentralRow-1;i>=0;i--) {
      xbin[i+i_skip]=xbin[i+1+i_skip]-celws[NXCellCentralRow-i-1];
   }
   i++;  // go back to the last item entered
   if (xmin < xbin[i]) {xbin[0] = xmin;nbin++;}
   int nbin_span=0;
   for(i=NXCellCentralRow-1;i>=0;i--) {
      if (fabs(xbin[i])>DiamondUpperCellW) break;
      nbin_span++;
   }
   
   TH2F* tofocc = new TH2F("tofocc","ToF occupancy;X(mm);Y(mm)",nbin,const_cast<double*>(xbin),(int)DiamondNYCell,
                                          -DiamondNYCell*DiamondCellHeight/2.,DiamondNYCell*DiamondCellHeight/2.);
   if (arm!="F"&&arm!="B") return;
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      cout << "Tree not found."<<endl;
      return;
   }
// new method
   int lat_1st_bin=tofocc->GetXaxis()->FindBin(ToFXPos-DiamondUpperCellW)+1;
   int lat_lst_bin=tofocc->GetXaxis()->FindBin(ToFXPos)-1;
//
   std::vector<double> tofX;
   std::vector<double> tofY;
   int nhits=0;
   TBranch* b_tofX=0;
   TBranch* b_tofY=0;
   TBranch* b_nhits=0;
 
   T->SetMakeClass(1);
   T->SetBranchAddress(Form("%s.Arm%s.ToFDet.X",fRecoBranchName.c_str(),arm.c_str()),&tofX,&b_tofX);
   T->SetBranchAddress(Form("%s.Arm%s.ToFDet.Y",fRecoBranchName.c_str(),arm.c_str()),&tofY,&b_tofY);
   T->SetBranchAddress(Form("%s.Arm%s.ToFDet.NHits",fRecoBranchName.c_str(),arm.c_str()),&nhits,&b_nhits);
   int nentries=T->GetEntries();
   for(i=0;i<nentries;i++) {
      Long64_t ientry = T->LoadTree(i);
      if (ientry<0) break;
      nhits=0;tofX.clear();tofY.clear();
      b_nhits->GetEntry(ientry);
      b_tofX->GetEntry(ientry);
      b_tofY->GetEntry(ientry);
      if (tofX.size()!=tofY.size()) {
         cout << "WARNING: screwed event" << endl;
         continue;
      }
      for(int j=0;j<nhits;j++) {
         int celId = ToFCellId(tofX.at(j),tofY.at(j));
         if (celId==0) continue;
         double xtof=ToFXFromCell(celId);
         double ytof=ToFYFromCell(celId);
         int biny = int(celId/100);
         int binx = celId-biny*100;
         if (biny==1||biny==3) {
/*
            for(int n=0;n<=nbin_span+1;n++) {
               tofocc->AddBinContent(tofocc->GetBin(NXCellCentralRow-n-1,biny));
            }
*/
            for(int n=lat_1st_bin;n<=lat_lst_bin;n++) tofocc->AddBinContent(tofocc->GetBin(n,biny));

	 } else {
           tofocc->Fill(xtof,ytof);
         }
      }
   }
   T->ResetBranchAddresses();
   T->SetMakeClass(0);
   tofocc->Scale(1./nentries);
}
void setToFParameters(std::string geo,double xpos)
{
     if (!fInitialized) setDataFormat();
     if (BeamX.count(1)==0||BeamX.count(2)==0||BeamX.count(3)==0||
        BeamY.count(1)==0||BeamY.count(2)==0||BeamY.count(3)==0||
        BeamX_RMS.count(1)==0||BeamX_RMS.count(2)==0||BeamX_RMS.count(3)==0||
        BeamY_RMS.count(1)==0||BeamY_RMS.count(2)==0||BeamY_RMS.count(3)==0) {
        std::cout << "ERROR: Beam parameters not given..." << std::endl;
        return;
     }
     ToFXPos=-(xpos*BeamX_RMS[3]+RPWindowThickness);
     setToFGeometry(geo);
}
void drawQuarticOccupancy(const string& arm,double xmin, double xmax, double ymin, double ymax)
{
   TH2F *tofocc = new TH2F("tofocc",Form("ToF %s @ %3.1fmm;X(mm);Y(mm)",ToFGeometryType.c_str(),ToFXPos),
                      5,ToFXPos-(NXCell)*QuarticCellWidth,ToFXPos,5,-NYCell*QuarticCellHeight/2.,(NYCell)*QuarticCellHeight/2.);
   if (arm!="F"&&arm!="B") return;
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      cout << "Tree not found."<<endl;
      return;
   }
   std::vector<double> tofX;
   std::vector<double> tofY;
   int nhits=0;
   TBranch* b_tofX=0;
   TBranch* b_tofY=0;
   TBranch* b_nhits=0;
 
   T->SetMakeClass(1);
   T->SetBranchAddress(Form("%s.Arm%s.ToFDet.X",fRecoBranchName.c_str(),arm.c_str()),&tofX,&b_tofX);
   T->SetBranchAddress(Form("%s.Arm%s.ToFDet.Y",fRecoBranchName.c_str(),arm.c_str()),&tofY,&b_tofY);
   T->SetBranchAddress(Form("%s.Arm%s.ToFDet.NHits",fRecoBranchName.c_str(),arm.c_str()),&nhits,&b_nhits);
   int nentries=T->GetEntries();
   for(int i=0;i<nentries;i++) {
      Long64_t ientry = T->LoadTree(i);
      if (ientry<0) break;
      nhits=0;tofX.clear();tofY.clear();
      b_nhits->GetEntry(ientry);
      b_tofX->GetEntry(ientry);
      b_tofY->GetEntry(ientry);
      if (tofX.size()!=tofY.size()) {
         cout << "WARNING: screwed event" << endl;
         continue;
      }
      for(int j=0;j<nhits;j++) {
         int celId = ToFCellId(tofX.at(j),tofY.at(j));
         double xtof=ToFXFromCell(celId);
         double ytof=ToFYFromCell(celId);
         tofocc->Fill(xtof,ytof);
      }
   }
   T->ResetBranchAddresses();
   T->SetMakeClass(0);
   tofocc->Scale(1./nentries);
}
void drawToFTotalMultiplicity(const std::string& arm)
{
   if (arm!="F"&&arm!="B") return;
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      cout << "Tree not found."<<endl;
      return;
   }
   TH1F *tofMultiplicity = new TH1F("tofMultiplicity","ToF total multiplicity/event;N hits",9,1,10);
   std::vector<double> tofX;
   std::vector<double> tofY;
   int nhits=0;
   TBranch* b_tofX=0;
   TBranch* b_tofY=0;
   TBranch* b_nhits=0;
 
   T->SetMakeClass(1);
   T->SetBranchAddress(Form("%s.Arm%s.ToFDet.X",fRecoBranchName.c_str(),arm.c_str()),&tofX,&b_tofX);
   T->SetBranchAddress(Form("%s.Arm%s.ToFDet.Y",fRecoBranchName.c_str(),arm.c_str()),&tofY,&b_tofY);
   T->SetBranchAddress(Form("%s.Arm%s.ToFDet.NHits",fRecoBranchName.c_str(),arm.c_str()),&nhits,&b_nhits);
   int nentries=T->GetEntries();
   for(int i=0;i<nentries;i++) {
      Long64_t ientry = T->LoadTree(i);
      if (ientry<0) break;
      b_nhits->GetEntry(ientry);
      b_tofX->GetEntry(ientry);
      b_tofY->GetEntry(ientry);
      int tof_multi=0;
      for(int j=0;j<nhits;j++) {
         int cellid = ToFCellId(tofX.at(j),tofY.at(j));
         if (cellid>0) tof_multi++;
      }
      tofMultiplicity->Fill(tof_multi);
   }
   T->ResetBranchAddresses();
   T->SetMakeClass(0);
   tofMultiplicity->Scale(1./nentries);
}
void drawToFCellMultiplicity(const std::string& arm)
{
   if (arm!="F"&&arm!="B") return;
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      cout << "Tree not found."<<endl;
      return;
   }
   TH1F *tofCellMultiplicity = new TH1F("tofCellMultiplicity","ToF Cell multiplicity",9,1,10);
   std::vector<double> tofX;
   std::vector<double> tofY;
   int nhits=0;
   TBranch* b_tofX=0;
   TBranch* b_tofY=0;
   TBranch* b_nhits=0;

   T->SetMakeClass(1);
   T->SetBranchAddress(Form("%s.Arm%s.ToFDet.X",fRecoBranchName.c_str(),arm.c_str()),&tofX,&b_tofX);
   T->SetBranchAddress(Form("%s.Arm%s.ToFDet.Y",fRecoBranchName.c_str(),arm.c_str()),&tofY,&b_tofY);
   T->SetBranchAddress(Form("%s.Arm%s.ToFDet.NHits",fRecoBranchName.c_str(),arm.c_str()),&nhits,&b_nhits);
   int nentries=T->GetEntries();
   for(int i=0;i<nentries;i++) {
      Long64_t ientry = T->LoadTree(i);
      if (ientry<0) break;
      b_nhits->GetEntry(ientry);
      b_tofX->GetEntry(ientry);
      b_tofY->GetEntry(ientry);
      std::map<int,int> cells;
      cells.clear();
      for(int j=0;j<nhits;j++) {
         int cellid = ToFCellId(tofX.at(j),tofY.at(j));
         if (cellid>0) cells[cellid]=1;
      }
      tofCellMultiplicity->Fill(cells.size());
   }
   T->ResetBranchAddresses();
   T->SetMakeClass(0);
   tofCellMultiplicity->Scale(1./nentries);
}
void drawToFHitPerCellMultiplicity(const std::string& arm)
{
}
void setBeamParameters(int pos, double x,double y,double rms_x,double rms_y)
{
     if (pos<1||pos>3) {
        cout << "ERROR: unknown detector position index: " << pos << endl;
        cout << "       It should be 1 for Tracker 1, 2 for Tracker 2 and 3 for ToF."<<endl;
        return;
     }
     BeamX[pos]=x;
     BeamY[pos]=y;
     BeamX_RMS[pos]=rms_x;
     BeamY_RMS[pos]=rms_y;
}
void  get_t_atMCframe(double cr_ang)
{
   double fBeamEnergy = 6500.;
   double ProtonMassSQ = pow(0.93827,2);
   double fBeamMomentum = sqrt(fBeamEnergy*fBeamEnergy - ProtonMassSQ);
   double microrad = 1.e-6;
   TTree* T = (TTree*)gDirectory->FindObject(fTreeName.c_str());
   if (!T) {
      cout << "Tree not found."<<endl;
      return;
   }
   std::vector<double> ArmF_t;
   std::vector<double> ArmF_xi;
   std::vector<double> ArmF_p;
   std::vector<double> ArmF_phi;
   std::vector<double> ArmF_theta;
   int nvtx=0;
   TBranch* b_ArmF_t=0;
   TBranch* b_ArmF_xi=0;
   TBranch* b_ArmF_p=0;
   TBranch* b_ArmF_phi=0;
   TBranch* b_ArmF_theta=0;
   //TBranch* b_nvtx=0;

   T->SetMakeClass(1);
   T->SetBranchAddress(Form("%s.ArmF.t",fRecoBranchName.c_str()),&ArmF_t,&b_ArmF_t);
   T->SetBranchAddress(Form("%s.ArmF.xi",fRecoBranchName.c_str()),&ArmF_xi,&b_ArmF_xi);
   T->SetBranchAddress(Form("%s.ArmF.momentum",fRecoBranchName.c_str()),&ArmF_p,&b_ArmF_p);
   T->SetBranchAddress(Form("%s.ArmF.phi",fRecoBranchName.c_str()),&ArmF_phi,&b_ArmF_phi);
   T->SetBranchAddress(Form("%s.ArmF.theta",fRecoBranchName.c_str()),&ArmF_theta,&b_ArmF_theta);
   //T->SetBranchAddress(Form("%s.ArmF.Nvtx",fRecoBranchName.c_str()),&nvtx,&b_nvtx);
   int nentries=T->GetEntries();
   for(int i=0;i<nentries;i++) {
      Long64_t ientry = T->LoadTree(i);
      if (ientry<0) break;
      //b_nvtx->GetEntry(ientry);
      b_ArmF_t->GetEntry(ientry);
      if (ArmF_t.at(0)==0) continue;
      b_ArmF_xi->GetEntry(ientry);
      b_ArmF_p->GetEntry(ientry);
      b_ArmF_phi->GetEntry(ientry);
      b_ArmF_theta->GetEntry(ientry);
      double px = ArmF_p.at(0)*(sin(ArmF_theta.at(0))*cos(ArmF_phi.at(0))-cos(cr_ang*microrad));
      double py = ArmF_p.at(0)*(sin(ArmF_theta.at(0))*sin(ArmF_phi.at(0)));
      double pz = ArmF_p.at(0)*(cos(ArmF_theta.at(0)-cr_ang*microrad));
      double mom = sqrt(px*px+py*py+pz*pz);
      double theta = atan(py/pz);
      double energy  = sqrt(px*px+py*py+pz*pz+ProtonMassSQ);
      double t  = -2.*(ProtonMassSQ-fBeamEnergy*energy+fBeamMomentum*mom*cos(theta));
      double xi = (1.-energy/fBeamEnergy);
      cout << "P: " << ArmF_p.at(0) << " theta: " << ArmF_theta.at(0) << " t: " << ArmF_t.at(0) << " xi: " << ArmF_xi.at(0) << endl;
      cout << mom << " " << theta << " " << t << " " << xi << endl;
      cout << "\n" << endl;
   }
}
