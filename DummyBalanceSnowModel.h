#ifndef DUMMY_BALANCE_SNOW_MODEL
#define DUMMY_BALANCE_SNOW_MODEL

#include <datastr/GeographicMap.h>

using namespace std;

namespace openworld {
  class DummyBalanceSnowModel : public SnowModel {
  protected:
    MatrixGeographicMap<double> volumes;

  public:
  DummyBalanceSnowModel(GeographicMap<double>& initial, DividedRange time)
    : SnowModel(time), volumes(initial) {
    }

    virtual SnowModel* clone() {
      return new DummyBalanceSnowModel(volumes, time);
    }

    // returns the snowcover
    virtual GeographicMap<double>& operator[](Measure tt) {
      return 100 * (volumes > 0);
    }

    virtual void inform(GeographicMap<double>& newMeltVolume, GeographicMap<double>& newSnowVolume) {
      volumes += newSnowVolume - newMeltVolume;
    }
  };
}

#endif
