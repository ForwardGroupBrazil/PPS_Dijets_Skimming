#include "DataFormats/PPSObjects/interface/PPSSpectrometer.h"
#include "DataFormats/PPSObjects/interface/PPSData.h"
#include "DataFormats/PPSObjects/interface/PPSDetector.h"
#include "DataFormats/Common/interface/Wrapper.h"

//namespace { struct dictionary {
    PPSData         PpsData;
    PPSDetector     PpsDetector;
    PPSSpectrometer PpsSpectrometer;
    edm::Wrapper<PPSSpectrometer> dummy;
//};}
