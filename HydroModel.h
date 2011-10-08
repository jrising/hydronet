#ifndef HYDROMODEL_H
#define HYDROMODEL_H

#include <datastr/GeographicMap.h>
#include <datastr/ScaledGeographicMap.h>
#include <datastr/TemporalGeographicMap.h>
#include <datastr/PartialConfidenceTemporalGeographicMap.h>
#include <datastr/TimeSeries.h>
#include <datastr/DividedRange.h>

using namespace openworld;

class HydroModel {
 protected:
  DividedRange latitudes;
  DividedRange longitudes;

  // Inputs
  GeographicMap<double>* elevation;
  GeographicMap<double>* slope;
  GeographicMap<double>* direction;

  PartialConfidenceTemporalGeographicMap<double>* precipitation;
  PartialConfidenceTemporalGeographicMap<double>* surfaceTemp;
  TemporalGeographicMap<double>* snowCover;

  double precipCoefficient;
  double meltCoefficient;

  // Results (evaluated at time = now)
  Measure now;
  MatrixGeographicMap<double>* riverVolume; // m^3
  MatrixGeographicMap<double>* surfaceVolume; // m^3
  MatrixGeographicMap<double>* riverConf;
  MatrixGeographicMap<double>* surfaceConf;  

  // Records of outflows from a particular cell
  unsigned rrFlow, ccFlow;
  list<double> outFlow;

  double minSurface, minRiver;
  double maxSurfaceVelocity, maxRiverVelocity;

 public:
  HydroModel(DividedRange latitudes, DividedRange longitudes, double precipCoefficient, double meltCoefficient, double latFlow, double longFlow);
  HydroModel(DividedRange latitudes, DividedRange longitudes, double precipCoefficient, double meltCoefficient, unsigned rrFlow, unsigned ccFlow);

  void setElevation(GeographicMap<double>* elevation);
  void setDInfinity(GeographicMap<double>* slope, GeographicMap<double>* direction);
  void setPrecipitation(PartialConfidenceTemporalGeographicMap<double>* precipitation);
  void setTemperature(PartialConfidenceTemporalGeographicMap<double>* surfaceTemp);
  void setSnowCover(TemporalGeographicMap<double>* snowCover);

  void runTo(time_t time);
  void stepDay();

  list<double> getOutFlows();
  GeographicMap<double>* getRiverVolume();

  // Utility functions

  double calcManning(double coeff, double radius, double slope);
  GeographicMap<double>& calcManningRiver(GeographicMap<double>* volume);
  GeographicMap<double>& calcManningSurface(GeographicMap<double>* volume);

  double calcMinimumDistance(GeographicMap<double>& map);
  double calcRiverWidth(double volume);
};

#endif
