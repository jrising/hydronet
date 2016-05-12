#ifndef COMPARE_SNOW_MODEL
#define COMPARE_SNOW_MODEL

#include "BackupSnowModel.h"
#include <datastr/ConstantGeographicMap.h>

using namespace std;

namespace openworld {
  class CompareSnowModel : public SnowModel {
  protected:
    BackupSnowModel& backup;
    DividedRange backuptime;
    int lastindex;

    SnowModel& compare;
    string outfile;

  public:
  CompareSnowModel(BackupSnowModel& backup, SnowModel& compare, DividedRange backuptime, string outfile)
    : SnowModel(backup.getTimes()), backup(backup), backuptime(backuptime), compare(compare) {
      this->lastindex = -1;
      this->outfile = outfile;
    }

    virtual SnowModel* clone() {
      CompareSnowModel* clone = new CompareSnowModel(*((BackupSnowModel*) backup.clone()), *compare.clone(), backuptime, outfile);
      clone->lastindex = lastindex;

      return clone;
    }

    // returns the snowcover
    virtual GeographicMap<double>& operator[](Measure tt) {
      GeographicMap<double>& result = backup[tt];

      int index = backuptime.inRange(tt);
      if (index >= 0 && index != lastindex) {
        GeographicMap<double>& comparison = compare[tt];

        cout << "SNOW COMPARE" << endl;
        fstream fs;
        fs.open(outfile.c_str(), fstream::out | fstream::app);
        for (unsigned rr = 0; rr < comparison.getLatitudes().count(); rr++)
          for (unsigned cc = 0; cc < comparison.getLongitudes().count(); cc++) {
            Measure latitude(Inds::lat), longitude(Inds::lon);
            comparison.calcLatitudeLongitude(rr, cc, latitude, longitude);
            if (result.validLocation(latitude, longitude))
              fs << rr << "\t" << cc << "\t" << result.getDouble(latitude, longitude) << "\t" << comparison.getCellConst(rr, cc) << endl;
          }
        fs.close();

        lastindex = backuptime.inRange(tt);
      }

      return result;
    }

    virtual void inform(GeographicMap<double>& newMeltVolume, GeographicMap<double>& newSnowVolume) {
      if (lastindex >= 0)
        compare.inform(newMeltVolume, newSnowVolume);
    }
  };
}

#endif
