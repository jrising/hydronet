#ifndef BALANCE_SNOW_MODEL
#define BALANCE_SNOW_MODEL

#include <datastr/GeographicMap.h>
#include <datastr/ConstantGeographicMap.h>

using namespace std;

namespace openworld {
  class BalanceSnowModel : public SnowModel {
  protected:
    MatrixGeographicMap<double> volumes;
    time_t initialTime;
    
    double lastCoverCoefficient;
    double lastSnowsCoefficient;
    double lastMeltsCoefficient;
    double lastMassCoefficient;

    int lastIndex;
    MatrixGeographicMap<double> sinceLastSnows;
    MatrixGeographicMap<double> sinceLastMelts;
    MatrixGeographicMap<double> lastCover;

  public:
  BalanceSnowModel(GeographicMap<double>& initialVolumes, GeographicMap<double>& initialCover, DividedRange time, double lastCoverCoefficient, double lastSnowsCoefficient, double lastMeltsCoefficient, double lastMassCoefficient)
    : SnowModel(time), volumes(initialVolumes), sinceLastSnows(ConstantGeographicMap<double>(initialVolumes.getLatitudes(), initialVolumes.getLongitudes(), 0)), sinceLastMelts(ConstantGeographicMap<double>(initialVolumes.getLatitudes(), initialVolumes.getLongitudes(), 0)), lastCover(initialCover) {
      this->lastCoverCoefficient = lastCoverCoefficient;
      this->lastSnowsCoefficient = lastSnowsCoefficient;
      this->lastMeltsCoefficient = lastMeltsCoefficient;
      this->lastMassCoefficient = lastMassCoefficient;

      lastIndex = 0;
    }

    virtual SnowModel* clone() {
      return new BalanceSnowModel(volumes, lastCover, time, lastCoverCoefficient, lastSnowsCoefficient,
                                  lastMeltsCoefficient, lastMassCoefficient);
    }

    virtual DividedRange getLongitudes() {
      return volumes.getLongitudes();
    }

    virtual DividedRange getLatitudes() {
      return volumes.getLatitudes();
    }

    // returns the snowcover
    virtual GeographicMap<double>& operator[](Measure tt) {
      int tIndex = time.inRange(tt);
      if (tIndex == -1)
        throw runtime_error("Time is before start in BalanceSnowModel");

      if (tIndex != lastIndex) {
        /*cout << "SNOW COMPARE" << endl;
        for (unsigned rr = 0; rr < volumes.getLatitudes().count(); rr++)
          for (unsigned cc = 0; cc < volumes.getLongitudes().count(); cc++) {
            cout << rr << ", " << cc << ": " << debugInfo(rr, cc) << endl;
            }*/

        lastCover = lastCoverCoefficient*lastCover + lastSnowsCoefficient*sinceLastSnows + lastMeltsCoefficient*sinceLastMelts + lastMassCoefficient*volumes;
        lastCover = lastCover * ((lastCover > 0) * (lastCover < 100)) + 100 * (lastCover > 100);

        volumes += sinceLastSnows - sinceLastMelts;
        volumes = volumes * (volumes > 0);
        sinceLastMelts.loadConstantInto(0);
        sinceLastSnows.loadConstantInto(0);

        lastIndex = tIndex;
      }

      return lastCover;
    }

    virtual string debugInfo(unsigned rr, unsigned cc) {
      stringstream out;
      out << lastCover.getCellConst(rr, cc) << "\t" << sinceLastSnows.getCellConst(rr, cc) << "\t" << sinceLastMelts.getCellConst(rr, cc) << "\t" << volumes.getCellConst(rr, cc);
      return out.str();
    }

    virtual void inform(GeographicMap<double>& newMeltVolume, GeographicMap<double>& newSnowVolume) {
      sinceLastMelts += newMeltVolume;
      sinceLastSnows += newSnowVolume;
    }
  };
}

#endif
