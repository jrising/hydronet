#ifndef SNOW_MODEL
#define SNOW_MODEL

#include <datastr/PartialConfidenceTemporalGeographicMap.h>
#include <datastr/TemporalGeographicMap.h>
#include <time.h>

using namespace std;

namespace openworld {
  class SnowModel : public PartialConfidenceTemporalGeographicMap<double> {
  public:
    SnowModel(DividedRange time)
      : PartialConfidenceTemporalGeographicMap<double>(time) {
    }

    virtual SnowModel* clone() = 0;
    virtual GeographicMap<double>& operator[](Measure tt) = 0;
    // needs to be on same grid as snowmodel
    virtual void inform(GeographicMap<double>& newMeltVolume, GeographicMap<double>& newSnowVolume) = 0;
    virtual string debugInfo(unsigned rr, unsigned cc) {
      return "";
    }

    /* Can't scale, because volumes need to be combined...
    virtual void scaledInform(GeographicMap<double>& newMeltVolume, GeographicMap<double>& newSnowVolume, GeographicMap<double>& scaleTo) {
      ScaledGeographicMap<double> scaledNewMeltVolume(newMeltVolume, scaleTo.getLatitudes(), scaleTo.getLongitudes(), 0.0);
      ScaledGeographicMap<double> scaledNewSnowVolume(newSnowVolume, scaleTo.getLatitudes(), scaleTo.getLongitudes(), 0.0);
      inform(scaledNewMeltVolume, scaledNewSnowVolume);
    }*/
  };
}

#endif
