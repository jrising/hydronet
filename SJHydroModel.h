#ifndef HYDROMODEL_H
#define HYDROMODEL_H

#include <datastr/TimeSeries.h>
#include <datastr/DividedRange.h>
#include "SnowModel.h"

#define ZERO_CELSIUS 273.15
#define RAIN_ALL_TEMPERATURE (ZERO_CELSIUS + 2)
#define SNOW_ALL_TEMPERATURE ZERO_CELSIUS

using namespace openworld;

template <class TBase, class TNumeric, class TLogical>
class SJHydroModel {
 protected:
  Indicator timeind;

  int verbose;
  
  double precipMult;
  double tempAdd;
  double snowDiff;

  // Results of internalStepDay
  // Updated sequentially
  TNumeric& currentGroundRainVolume;
  TNumeric& currentGroundMeltVolume;
  TNumeric& currentGroundConf;

  // Refreshed every step
  TNumeric& newRainVolume;
  TNumeric& newMeltVolume;
  TNumeric& newVolumeConf;
  
  list<double> outFlowRain;
  list<double> outFlowMelt;
  Measure now;

 public:
  SJHydroModel(Indicator timeind, double meltDegreeDayFactor, double meltDegreeDaySlope, double rainRunoffCoefficient, double meltRunoffCoefficient, double groundCoefficient, double groundToBaseflowDay, double rainOnSnowCoefficient, double surfaceEvaporationFactor, double riverEvaporationFactor, TNumeric& initialGroundRainVolume, TNumeric& initialGroundMeltVolume, TNumeric& initialGroundConf);
  SJHydroModel(SJHydroModel& copy);

  virtual ~SJHydroModel();

  double meltDegreeDayFactor;
  double meltDegreeDaySlope;
  double rainRunoffCoefficient;
  double meltRunoffCoefficient;
  double groundCoefficient;
  double groundToBaseflowDay;
  double rainOnSnowCoefficient;
  double surfaceEvaporationFactor;
  double riverEvaporationFactor;

  void setPrecipitationScaling(double precipMult = 1);
  void setTemperatureAddition(double tempAdd = 0);
  void setSnowCoverDifference(double snowDiff = 0);

  virtual time_t getTime();
  virtual void runTo(time_t time);
  virtual void stepDay() = NULL;
  virtual void internalStepDay() = NULL;
  
  list<double> getOutFlowsRain();
  list<double> getOutFlowsMelt();
  unsigned getOutFlowsCount();

  // Serializable protocol

  friend istream& operator>>(istream& in, SJHydroModel& sink) {
    sink.timeind = Indicator::streamExtract(in);
    in >> sink.precipMult >> sink.tempAdd >> sink.snowDiff;

    // DO NOT INPUT CURRENT STATE
    return in;
  }

  friend ostream& operator<<(ostream& os, SJHydroModel& source) {
    source.timeind.streamInsert(os);
    os << source.precipMult << " " << source.tempAdd << " " << source.snowDiff;

    // DO NOT OUTPUT CURRENT STATE
    return os;
  }

 protected:
  TNumeric& weightedFraction(TNumeric& weights1, TNumeric& values1, TNumeric& weights2, TNumeric& values2);
};

#endif
