#ifndef PPSSpectrometer_h
#define PPSSpectrometer_h
#include <vector>
#include <utility>
#include "TObject.h"
#include "DataFormats/PPSObjects/interface/PPSData.h"

class PPSSpectrometer: public TObject {
public:
      PPSSpectrometer();
      virtual ~PPSSpectrometer(){};

      int            Nvtx;
      std::vector<double> vtxX;
      std::vector<double> vtxY;
      std::vector<double> vtxZ;
      std::vector<std::pair<unsigned int,unsigned int> > vtxTracks;
      PPSData ArmF;
      PPSData ArmB;

      int  AddVertex(double vtx,double vty, double vtz) {
                    vtxX.push_back(vtx);
                    vtxY.push_back(vty);
                    vtxZ.push_back(vtz);
                    Nvtx++;
                    return Nvtx;
           };

      void clear() { Nvtx=0;vtxX.clear(); vtxY.clear(); vtxZ.clear();vtxTracks.clear();ArmF.clear(); ArmB.clear(); };

ClassDef(PPSSpectrometer,1);
};
#endif
