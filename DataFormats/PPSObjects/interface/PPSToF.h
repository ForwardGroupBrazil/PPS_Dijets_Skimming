#ifndef PPSToF_h
#define PPSToF_h
#include <vector>
#include "TObject.h"

class PPSToF: public TObject {
      public:
            PPSToF();
            virtual ~PPSToF() {};

      public:
            const int                 NCellX;
            const int                 NCellY;
            const double              DetW; // width (X, horizontal dimension in mm)
            const double              DetH; // height(Y, vertical dimension in mm)
//
            int                 DetId;
            int                 NHits;
            std::vector<double> CelId; // Starts numbering from 1 in the bottom (-Y) closest corner to the beam
            std::vector<double> ToF;
            void clear() {DetId=0;NHits=0;CelId.clear();ToF.clear();};
            void AddHit(double x, double y, double z, double thx = 0., double thy=0.,int stop=0)
                 {
                 };
ClassDef(PPSToF,1);
};
#endif
