#ifndef HYDRONETMODEL_H
#define HYDRONETMODEL_H

#include <datastr/GeographicMap.h>
#include <datastr/ScaledGeographicMap.h>
#include <datastr/TemporalGeographicMap.h>
#include <datastr/PartialConfidenceTemporalGeographicMap.h>
#include <datastr/TimeSeries.h>
#include <datastr/DividedRange.h>
#include "DInfinityMap.h"

using namespace openworld;

class HydroNetModel {
 protected:
  DividedRange latitudes;
  DividedRange longitudes;

  // Inputs
  HydroNet* net;
  HydroOutputNode* out;

  PartialConfidenceTemporalGeographicMap<double>* precipitation;
  PartialConfidenceTemporalGeographicMap<double>* surfaceTemp;
  TemporalGeographicMap<double>* snowCover;

  double precipCoefficient;
  double meltCoefficient;

  // Results (evaluated at time = now)
  Measure now;

 public:
  HydroNetModel(DividedRange latitudes, DividedRange longitudes, double precipCoefficient, double meltCoefficient);

  void setup(GeographicMap<bool>* mask_coarse, GeographicMap<double>* mask_fine, GeographicMap<double>* slope_fine, DInfinityMap* direction_fine);

  void setPrecipitation(PartialConfidenceTemporalGeographicMap<double>* precipitation);
  void setTemperature(PartialConfidenceTemporalGeographicMap<double>* surfaceTemp);
  void setSnowCover(TemporalGeographicMap<double>* snowCover);

  void runTo(time_t time);
  void stepDay();

  list<double> getOutFlows();
};

#endif
