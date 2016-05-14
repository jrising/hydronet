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

class SJHydroNetModel;

class SJHydroNetModelStepCallback {
 public:
  virtual void setup(SJHydroNetModel& model) = 0;
  virtual void post(SJHydroNetModel& model) = 0;
};

class SJHydroNetModel {
 protected:
  DividedRange latitudes;
  DividedRange longitudes;
  Indicator timeind;
  int verbose;

  // Inputs
  HydroNet* net;
  HydroOutputNode* out;

  PartialConfidenceTemporalGeographicMap<double>* precipitation;
  PartialConfidenceTemporalGeographicMap<double>* surfaceTemp;
  SnowModel* snowModel;
  GeographicMap<double>* elevation;
  GeographicMap<double>* mmdayToVolume;
  GeographicMap<double>* mask_coarse;

  SJHydroNetModelStepCallback* stepCallback;

  double precipMult;
  double tempAdd;
  double snowDiff;

  // Results (evaluated at time = now)
  GeographicMap<double>* currentGroundRainVolume;
  GeographicMap<double>* currentGroundMeltVolume;
  GeographicMap<double>* currentGroundConf;

  list<double> outFlowRain;
  list<double> outFlowMelt;
  Measure now;

 public:
  SJHydroNetModel(DividedRange latitudes, DividedRange longitudes, Indicator timeind, double meltDegreeDayFactor, double meltDegreeDaySlope, double rainRunoffCoefficient, double meltRunoffCoefficient, double groundCoefficient, double groundToBaseflowDay, double rainOnSnowCoefficient, double surfaceEvaporationFactor, double riverEvaporationFactor);
  SJHydroNetModel(SJHydroNetModel& copy);

  ~SJHydroNetModel();

  double meltDegreeDayFactor;
  double meltDegreeDaySlope;
  double rainRunoffCoefficient;
  double meltRunoffCoefficient;
  double groundCoefficient;
  double groundToBaseflowDay;
  double rainOnSnowCoefficient;
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
  void setVerbosity(int verbose);

  time_t getTime();
  void runTo(time_t time);
  void stepDay();

  template <class TNumeric, class TLogical>
    static void stepDayHeight(TNumeric& scaledPrecipitation, TNumeric& scaledSurfaceTemp,
                              TNumeric& fracSnowCover, TNumeric& fullMeltDegreeDayFactor,
                              double rainRunoffCoefficient, double meltRunoffCoefficient,
                              double groundCoefficient, double rainOnSnowCoefficient,
                              TNumeric*& newSnowMeltHeightPtr, TNumeric*& newSnowAccumHeightPtr,
                              TNumeric*& newRainRunoffHeightPtr, TNumeric*& newMeltRunoffHeightPtr,
                              TNumeric*& newRainGroundHeightPtr, TNumeric*& newMeltGroundHeightPtr,
                              TNumeric*& newDirectHeightConfPtr, bool verbose) {
    // Add precipitation
    if (verbose)
      cout << "Determining valid inputs" << endl;
    TLogical& validPrecipitation = scaledPrecipitation >= 0.0;
    TLogical& validSurfaceTemp = scaledSurfaceTemp >= 100.0;
    TLogical& validSnowCover = fracSnowCover >= 0.0;
    TLogical& aboveZeroCelsius = scaledSurfaceTemp >= ZERO_CELSIUS;

    if (verbose)
      cout << "Adding surface flow" << endl;
    // rain-related calculations
    TNumeric& rainPortion = ((scaledSurfaceTemp - SNOW_ALL_TEMPERATURE) / (RAIN_ALL_TEMPERATURE - SNOW_ALL_TEMPERATURE)) * (scaledSurfaceTemp > SNOW_ALL_TEMPERATURE) * (scaledSurfaceTemp <= RAIN_ALL_TEMPERATURE) + (scaledSurfaceTemp > RAIN_ALL_TEMPERATURE);
    TNumeric& newRainAllHeight = rainPortion * scaledPrecipitation * (validPrecipitation * validSurfaceTemp);
    TNumeric& newSnowfreeRainHeight = newRainAllHeight * (1.0 - fracSnowCover) * validSnowCover;

    // snow-related calculations
    TNumeric& newRainOnSnowRainHeight = newRainAllHeight * fracSnowCover * validSnowCover; // uses melt runoff, because rain on snow like melt
    TNumeric& newRainOnSnowMeltHeight = rainOnSnowCoefficient * (scaledSurfaceTemp - ZERO_CELSIUS) * aboveZeroCelsius * newRainOnSnowRainHeight;
    //TNumeric& validMeltDegreeDayFactor = fullMeltDegreeDayFactor * (fullMeltDegreeDayFactor > 0);
    TNumeric& newDegreeDayMeltHeight = fullMeltDegreeDayFactor * fracSnowCover * (scaledSurfaceTemp - ZERO_CELSIUS) * aboveZeroCelsius * validSnowCover;

    TNumeric& newSnowMeltHeight = newDegreeDayMeltHeight + newRainOnSnowMeltHeight;
    TNumeric& newSnowAccumHeight = (1.0 - rainPortion) * scaledPrecipitation * (validPrecipitation * validSurfaceTemp);

    // final rain-related and melt-related heights
    TNumeric& newRainRunoffHeight = rainRunoffCoefficient * newSnowfreeRainHeight + meltRunoffCoefficient * newRainOnSnowRainHeight;
    TNumeric& newMeltRunoffHeight = meltRunoffCoefficient * newSnowMeltHeight;

    TNumeric& newRainGroundHeight = groundCoefficient * ((1 - rainRunoffCoefficient) * newSnowfreeRainHeight + (1 - meltRunoffCoefficient) * newRainOnSnowRainHeight);
    TNumeric& newMeltGroundHeight = groundCoefficient * (1 - meltRunoffCoefficient) * newSnowMeltHeight; // CHANGE from paper: don't use M, which is multiplied by melt coefficient

    TNumeric& newDirectHeightConf = *(tew_(TNumeric(validPrecipitation * validSurfaceTemp * validSnowCover)));

    newSnowMeltHeightPtr = &newSnowMeltHeight;
    newSnowAccumHeightPtr = &newSnowAccumHeight;
    newRainRunoffHeightPtr = &newRainRunoffHeight;
    newMeltRunoffHeightPtr = &newMeltRunoffHeight;
    newRainGroundHeightPtr = &newRainGroundHeight;
    newMeltGroundHeightPtr = &newMeltGroundHeight;
    newDirectHeightConfPtr = &newDirectHeightConf;
  }

  template <class TTemporal, class TNumeric, class TLogical>
    static void runToHeight(TTemporal& precipitation, TTemporal& surfaceTemp,
                            TTemporal& fracSnowCover, TNumeric& elevation,
                            double meltDegreeDayFactor, double meltDegreeDaySlope,
                            double rainRunoffCoefficient, double meltRunoffCoefficient,
                            double groundCoefficient, double rainOnSnowCoefficient,
                            TTemporal& newSnowMeltHeight, TTemporal& newSnowAccumHeight,
                            TTemporal& newRainRunoffHeight, TTemporal& newMeltRunoffHeight,
                            TTemporal& newRainGroundHeight, TTemporal& newMeltGroundHeight,
                            TTemporal& newDirectHeightConf,
                            Measure now, Measure endtime, bool verbose) {
  TNumeric& fullMeltDegreeDayFactor = (meltDegreeDayFactor + meltDegreeDaySlope * elevation);

  while (now.getValue() == 0 || now < endtime) {
    time_t now_time = now.getValue();
    struct tm* ptm = gmtime(&now_time);

    if (now.getValue() == 0) {
      // What is the earliest time we can handle?
      cout << "Start time: " << precipitation.getTimes().getMin() << ", " << surfaceTemp.getTimes().getMin() << ", " << fracSnowCover.getTimes().getMin() << endl;
      now = max(max(precipitation.getTimes().getMin(), surfaceTemp.getTimes().getMin()), fracSnowCover.getTimes().getMin());
      cout << "Starting at " << now << endl;
    } else {
      cout << ptm->tm_mday << "/" << ptm->tm_mon+1 << "/" << ptm->tm_year << ": " << now << " < " << endtime << endl;

      if (now.getIndicator() == Inds::unixtime)
        now += DividedRange::toTimespan(1);
      else
        now += 1/360.0;
    }

    TNumeric& nowPrecipitation = precipitation[now];
    TNumeric& nowSurfaceTemp = surfaceTemp[now];
    TNumeric& nowFracSnowCover = fracSnowCover[now];

    TNumeric *newSnowMeltHeightPtr, *newSnowAccumHeightPtr,
      *newRainRunoffHeightPtr, *newMeltRunoffHeightPtr,
      *newRainGroundHeightPtr, *newMeltGroundHeightPtr,
      *newDirectHeightConfPtr;
    stepDayHeight<TNumeric, TLogical>(nowPrecipitation, nowSurfaceTemp,
                                      nowFracSnowCover, fullMeltDegreeDayFactor,
                                      rainRunoffCoefficient, meltRunoffCoefficient,
                                      groundCoefficient, rainOnSnowCoefficient,
                                      newSnowMeltHeightPtr, newSnowAccumHeightPtr,
                                      newRainRunoffHeightPtr, newMeltRunoffHeightPtr,
                                      newRainGroundHeightPtr, newMeltGroundHeightPtr,
                                      newDirectHeightConfPtr, verbose);

    newSnowMeltHeight.add(now, *newSnowMeltHeightPtr);
    newSnowAccumHeight.add(now, *newSnowAccumHeightPtr);
    newRainRunoffHeight.add(now, *newRainRunoffHeightPtr);
    newMeltRunoffHeight.add(now, *newMeltRunoffHeightPtr);
    newRainGroundHeight.add(now, *newRainGroundHeightPtr);
    newMeltGroundHeight.add(now, *newMeltGroundHeightPtr);
    newDirectHeightConf.add(now, *newDirectHeightConfPtr);
  }
  }

  list<double> getOutFlowsRain();
  list<double> getOutFlowsMelt();
  unsigned getOutFlowsCount();

  HydroNet& getHydroNet() {
    return *net;
  }

  DividedRange getLatitudes() {
    return latitudes;
  }

  DividedRange getLongitudes() {
    return longitudes;
  }

  PartialConfidenceTemporalGeographicMap<double>& getPrecipitation() {
    return *precipitation;
  }

  PartialConfidenceTemporalGeographicMap<double>& getTemperature() {
    return *surfaceTemp;
  }

  SnowModel& getSnowModel() {
    return *snowModel;
  }

  GeographicMap<double>& getElevation() {
    return *elevation;
  }

  void setStepCallback(SJHydroNetModelStepCallback* stepCallback) {
    this->stepCallback = stepCallback;
  }

  SJHydroNetModelStepCallback* getStepCallback() {
    return stepCallback;
  }

  void resetStepCallback() {
    if (stepCallback)
      stepCallback->setup(*this);
  }

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
