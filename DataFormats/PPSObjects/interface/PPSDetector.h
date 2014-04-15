#ifndef PPSDetector_h
#define PPSDetector_h
#include <vector>
#include "TObject.h"

class PPSDetector: public TObject {
      public:
            PPSDetector();
            virtual ~PPSDetector() {};

      public:
            int                 DetId;
            int                 NHits;
            std::vector<double> X;
            std::vector<double> Y;
            std::vector<double> Z;
            std::vector<double> thetaX;
            std::vector<double> thetaY;
            std::vector<double> Eloss;
            std::vector<int>    HasStopped; // 1 if particle has stopped before Z position
            void clear() {DetId=0;NHits=0;X.clear();Y.clear();Z.clear();thetaX.clear();thetaY.clear();HasStopped.clear();};
            void AddHit(double x, double y, double z, double thx = 0., double thy=0.,int stop=0)
                 {X.push_back(x);Y.push_back(y);Z.push_back(z); HasStopped.push_back(stop);
                  thetaX.push_back(thx);thetaY.push_back(thy);NHits++;};
ClassDef(PPSDetector,1);
};
#endif
