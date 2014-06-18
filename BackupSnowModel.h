#ifndef BACKUP_SNOW_MODEL
#define BACKUP_SNOW_MODEL

#include <datastr/PartialConfidenceTemporalGeographicMap.h>
#include <datastr/TemporalGeographicMap.h>
#include <time.h>

using namespace std;

namespace openworld {
  class BackupSnowModel : public SnowModel {
  protected:
    TemporalGeographicMap<double>& data; // saves pointer
    TemporalGeographicMap<double>& backup; // saves pointer
    double backupConfidence;
    Measure backupStartTime;

  public:
    BackupSnowModel(TemporalGeographicMap<double>* data, TemporalGeographicMap<double>* backup,
                    DividedRange time, double backupConfidence, Measure backupStartTime)
      : SnowModel(time), data(*data), backup(*backup), backupStartTime(backupStartTime) {
      this->backupConfidence = backupConfidence;

      confs = new double[time.count()];
      for (unsigned ii = 0; ii < time.count(); ii++)
        confs[ii] = 0;
    }

    virtual SnowModel* clone() {
      return new BackupSnowModel(data.clone(), backup.clone(), this->time, backupConfidence, backupStartTime);
    }

    virtual GeographicMap<double>& operator[](Measure tt) {
      time_t timet = tt.getValue();
      struct tm tmstr = *gmtime(&timet);
      GeographicMap<double>& backupMap = backup[backupStartTime + DividedRange::toTimespan(tmstr.tm_yday)];

      if (tt >= data.getTimes().getMin() && tt < data.getTimes().getMax()) {
        GeographicMap<double>& dataMap = data[tt];

        GeographicMap<double>* result = tew_(MatrixGeographicMap<double>(dataMap.getLatitudes(), dataMap.getLongitudes()));
        unsigned valids = 0;
        for (unsigned rr = 0; rr < dataMap.getLatitudes().count(); rr++)
          for (unsigned cc = 0; cc < dataMap.getLongitudes().count(); cc++) {
            double val = dataMap.getCellConst(rr, cc);
            if (val >= 0) {
              valids++;
              result->getCell(rr, cc) = val;
            } else
              result->getCell(rr, cc) = backupMap.getCellConst(rr, cc);
          }

        confs[time.inRange(tt)] = backupConfidence + ((1 - backupConfidence) * valids) / (dataMap.getLatitudes().count() * dataMap.getLongitudes().count());
        return *result;
      } else {
        confs[time.inRange(tt)] = backupConfidence;
        return backupMap;
      }
    }

    virtual void inform(GeographicMap<double>& newMeltVolume, GeographicMap<double>& newSnowVolume) {
      // do nothing
    }
  };
}

#endif
