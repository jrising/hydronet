#ifndef DUMMY_COMPARE_SNOW_MODEL
#define DUMMY_COMPARE_SNOW_MODEL

#include "BackupSnowModel.h"
#include "DummyBalanceSnowModel.h"
#include <datastr/ConstantGeographicMap.h>

using namespace std;

namespace openworld {
  class DummyCompareSnowModel : public DummyBalanceSnowModel {
  protected:
    BackupSnowModel& backup;
    DividedRange backuptime;
    int lastindex;

  public:
  DummyCompareSnowModel(BackupSnowModel& backup, GeographicMap<double>& initial, DividedRange backuptime)
    : DummyBalanceSnowModel(initial, backup.getTimes()), backup(backup), backuptime(backuptime) {
      this->lastindex = -1;
    }

    virtual SnowModel* clone() {
      DummyCompareSnowModel* clone = new DummyCompareSnowModel(*((BackupSnowModel*) backup.clone()), volumes, backuptime);
      clone->lastindex = lastindex;

      return clone;
    }

    // returns the snowcover
    virtual GeographicMap<double>& operator[](Measure tt) {
      GeographicMap<double>& result = backup[tt];

      /*
      int index = backuptime.inRange(tt);
      if (index >= 0 && index != lastindex) {
        cout << "SNOW COMPARE," << tt << "," << index << endl;
        for (unsigned rr = 0; rr < volumes.getLatitudes().count(); rr++)
          for (unsigned cc = 0; cc < volumes.getLongitudes().count(); cc++) {
            Measure latitude(Inds::lat), longitude(Inds::lon);
            volumes.calcLatitudeLongitude(rr, cc, latitude, longitude);
            if (result.validLocation(latitude, longitude))
              cout << rr << "\t" << cc << "\t" << result.getDouble(latitude, longitude) << "\t" << volumes.getCellConst(rr, cc) << endl;
          }
        
        lastindex = backuptime.inRange(tt);
      }
      */

      return result;
    }
  };
}

#endif
