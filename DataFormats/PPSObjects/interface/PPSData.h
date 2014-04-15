#ifndef PPSData_h
#define PPSData_h
#include <vector>
#include "DataFormats/PPSObjects/interface/PPSDetector.h"

class PPSData {
  public:

      PPSData();
      virtual ~PPSData() {};

      int                 NTracks;
      std::vector<int>    vertexId;
      std::vector<double> momentum;
      std::vector<double> pT;
      std::vector<double> t;
      std::vector<double> xi;
      std::vector<double> eta;
      std::vector<double> phi;
      std::vector<double> theta;
      std::vector<double> XatTCL4;
      std::vector<double> XatTCL5;
      std::vector<double> ToF;
      std::vector<int>    stationId;
      std::vector<double> RatIP;
      PPSDetector TrkDet1;
      PPSDetector TrkDet2;
      PPSDetector ToFDet;
      PPSDetector ChkPoint;

      void Fill(int vtxid=-1,double _t=0,double _xi=0,double _p=0,double _pt=0,double _eta=0, double _phi=0,double _th=0,double _ImpPar=0,int st_id=0) {
               NTracks++; momentum.push_back(_p);pT.push_back(_pt);t.push_back(_t);xi.push_back(_xi);eta.push_back(_eta);
               phi.push_back(_phi);theta.push_back(_th);RatIP.push_back(_ImpPar);stationId.push_back(st_id);
           };
      void clear(){NTracks=0;
               momentum.clear(); t.clear(); xi.clear(); eta.clear(); phi.clear(); theta.clear(); 
               pT.clear(); ToF.clear(); stationId.clear();XatTCL4.clear();XatTCL5.clear();
               TrkDet1.clear();TrkDet2.clear();ToFDet.clear();ChkPoint.clear();RatIP.clear();
      };
ClassDef(PPSData,1);
};
#endif
