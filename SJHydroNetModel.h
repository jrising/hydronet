#ifndef HYDRONETMODEL_H
#define HYDRONETMODEL_H

#include <datastr/GeographicMap.h>
#include <datastr/ScaledGeographicMap.h>
#include <datastr/TemporalGeographicMap.h>
#include <datastr/PartialConfidenceTemporalGeographicMap.h>
#include <datastr/TimeSeries.h>
#include <datastr/DividedRange.h>
#include <tools/hydro/DInfinityMap.h>
#include "HydroNet.h"
#include "SnowModel.h"

#define ZERO_CELSIUS 273.15
#define RAIN_ALL_TEMPERATURE (ZERO_CELSIUS + 2)
#define SNOW_ALL_TEMPERATURE ZERO_CELSIUS

using namespace openworld;

class SJHydroNetModel {
 protected:
  DividedRange latitudes;
  DividedRange longitudes;
  Indicator timeind;

  // Inputs
  HydroNet* net;
  HydroOutputNode* out;

  PartialConfidenceTemporalGeographicMap<double>* precipitation;
  PartialConfidenceTemporalGeographicMap<double>* surfaceTemp;
  SnowModel* snowModel;
  GeographicMap<double>* elevation;
  GeographicMap<double>* mmdayToVolume;
  GeographicMap<double>* mask_coarse;

  double precipMult;
  double tempAdd;
  double snowDiff;

  // Results (evaluated at time = now)
  GeographicMap<double>* currentGroundRainVolume;
  GeographicMap<double>* currentGroundMeltVolume;
  GeographicMap<double>* currentGroundConf;

  string allcells;
  list<double> outFlowRain;
  list<double> outFlowMelt;
  Measure now;

 public:
  SJHydroNetModel(DividedRange latitudes, DividedRange longitudes, Indicator timeind, double meltDegreeDayFactor, double meltDegreeDaySlope, double rainRunoffCoefficient, double meltRunoffCoefficient, double groundCoefficient, double groundToBaseflowDay, double surfaceEvaporationFactor, double riverEvaporationFactor, string allcells = "");
  SJHydroNetModel(SJHydroNetModel& copy);

  ~SJHydroNetModel();

  double meltDegreeDayFactor;
  double meltDegreeDaySlope;
  double rainRunoffCoefficient;
  double meltRunoffCoefficient;
  double groundCoefficient;
  double groundToBaseflowDay;
  double surfaceEvaporationFactor;
  double riverEvaporationFactor;

  void setup(GeographicMap<double>* mask_coarse, GeographicMap<bool>* mask_fine, GeographicMap<double>* slope_fine, DInfinityMap* direction_fine, double mindist);

  void setPrecipitation(PartialConfidenceTemporalGeographicMap<double>* precipitation);
  void setTemperature(PartialConfidenceTemporalGeographicMap<double>* surfaceTemp);
  void setSnowModel(SnowModel* snowModel);
  void setElevation(GeographicMap<double>* elevation);

  void setPrecipitationScaling(double precipMult = 1);
  void setTemperatureAddition(double tempAdd = 0);
  void setSnowCoverDifference(double snowDiff = 0);

  time_t getTime();
  void runTo(time_t time);
  void stepDay();

  list<double> getOutFlowsRain();
  list<double> getOutFlowsMelt();
  unsigned getOutFlowsCount();

  // diagnostics
  list<pair<pair<pair<Measure, Measure>, pair<Measure, Measure> >, pair<bool, double> > > getAllEdges();

  // Serializable protocol

  friend istream& operator>>(istream& in, SJHydroNetModel& sink) {
    PointerReference reference;

    sink.latitudes = DividedRange::streamExtract(in);
    sink.longitudes = DividedRange::streamExtract(in);
    sink.timeind = Indicator::streamExtract(in);

    sink.net->streamExtract(in, reference);
    sink.out = HydroOutputNode::streamExtractPointer(in, reference);

    // DO NOT INPUT INPUT DATA

    in >> sink.precipMult >> sink.tempAdd >> sink.snowDiff;

    // DO NOT INPUT CURRENT STATE
    return in;
  }

  friend ostream& operator<<(ostream& os, SJHydroNetModel& source) {
      PointerTracker tracker;

      source.latitudes.streamInsert(os);
      source.longitudes.streamInsert(os);
      source.timeind.streamInsert(os);

      source.net->streamInsert(os, tracker);
      source.out->streamInsertPointer(os, tracker);

      // DO NOT OUTPUT INPUT DATA

      os << source.precipMult << " " << source.tempAdd << " " << source.snowDiff;

      // DO NOT OUTPUT CURRENT STATE
      return os;
  }
  //  friend istream& operator>>(istream& in, const SJHydroNetModel& sink);
  //  friend ostream& operator<<(ostream& os, const SJHydroNetModel& source);
  
 protected:
  GeographicMap<double>& weightedFraction(GeographicMap<double>& weights1, GeographicMap<double>& values1, GeographicMap<double>& weights2, GeographicMap<double>& values2);
};

#endif
