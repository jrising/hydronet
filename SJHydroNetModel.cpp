#include <algorithm>
#include "SJHydroNetModel.h"
#include <measure/Inds.h>
#include <memory/Transients.h>
#include <utils/Timer.h>
#include <datastr/ConstantGeographicMap.h>

SJHydroNetModel::SJHydroNetModel(DividedRange latitudes, DividedRange longitudes, Indicator timeind, double meltDegreeDayFactor, double meltDegreeDaySlope, double rainRunoffCoefficient, double meltRunoffCoefficient, double groundCoefficient, double groundToBaseflowDay, double rainOnSnowCoefficient, double surfaceEvaporationFactor, double riverEvaporationFactor)
  : latitudes(latitudes), longitudes(longitudes), timeind(timeind), now(0, timeind) {
  precipitation = NULL;
  surfaceTemp = NULL;
  snowModel = NULL;
  elevation = NULL;
  mask_coarse = NULL;

  stepCallback = NULL;

  precipMult = 1;
  tempAdd = 0;
  snowDiff = 0;

  now = 0;
  verbose = 1;

  this->meltDegreeDayFactor = meltDegreeDayFactor;
  this->meltDegreeDaySlope = meltDegreeDaySlope;
  this->rainRunoffCoefficient = rainRunoffCoefficient;
  this->meltRunoffCoefficient = meltRunoffCoefficient;
  this->groundCoefficient = groundCoefficient;
  this->groundToBaseflowDay = groundToBaseflowDay;
  this->rainOnSnowCoefficient = rainOnSnowCoefficient;
  this->surfaceEvaporationFactor = surfaceEvaporationFactor;
  this->riverEvaporationFactor = riverEvaporationFactor;

  currentGroundRainVolume = new MatrixGeographicMap<double>(latitudes, longitudes);
  currentGroundMeltVolume = new MatrixGeographicMap<double>(latitudes, longitudes);
  currentGroundConf = new MatrixGeographicMap<double>(latitudes, longitudes);

  // mm/day -> m^3/day
  ConstantGeographicMap<double> grid(latitudes, longitudes, 1);
  mmdayToVolume = &(grid.calcAreas() * .001);
  Transients::abandon(mmdayToVolume);
}

SJHydroNetModel::SJHydroNetModel(SJHydroNetModel& copy)
  : latitudes(copy.latitudes), longitudes(copy.longitudes), timeind(copy.timeind), now(copy.now) {

  map<HydroNode*, HydroNode*> translate;
  net = new HydroNet(*copy.net, translate);
  list<HydroNode*> ignore;
  out = (HydroOutputNode*) HydroNode::getCopy(copy.out, translate, ignore);

  precipitation = copy.precipitation->clone();
  surfaceTemp = copy.surfaceTemp->clone();
  snowModel = copy.snowModel->clone();
  elevation = copy.elevation; // don't delete!
  mask_coarse = copy.mask_coarse;

  precipMult = copy.precipMult;
  tempAdd = copy.tempAdd;
  snowDiff = copy.snowDiff;

  // direct copies
  mmdayToVolume = copy.mmdayToVolume;
  stepCallback = copy.stepCallback;

  meltDegreeDayFactor = copy.meltDegreeDayFactor;
  meltDegreeDaySlope = copy.meltDegreeDaySlope;
  rainRunoffCoefficient = copy.rainRunoffCoefficient;
  meltRunoffCoefficient = copy.meltRunoffCoefficient;
  groundCoefficient = copy.groundCoefficient;
  groundToBaseflowDay = copy.groundToBaseflowDay;
  rainOnSnowCoefficient = copy.rainOnSnowCoefficient;
  surfaceEvaporationFactor = copy.surfaceEvaporationFactor;
  riverEvaporationFactor = copy.riverEvaporationFactor;

  currentGroundRainVolume = new MatrixGeographicMap<double>(*copy.currentGroundRainVolume);
  currentGroundMeltVolume = new MatrixGeographicMap<double>(*copy.currentGroundMeltVolume);
  currentGroundConf = new MatrixGeographicMap<double>(*copy.currentGroundConf);

  outFlowRain = copy.outFlowRain;
  outFlowMelt = copy.outFlowMelt;
}

SJHydroNetModel::~SJHydroNetModel() {
  delete precipitation;
  delete surfaceTemp;
  delete snowModel;
  // don't delete elevation-- don't own it

  delete net;
  delete currentGroundRainVolume;
  delete currentGroundMeltVolume;
  delete currentGroundConf;
}

void SJHydroNetModel::setup(GeographicMap<double>* mask_coarse, GeographicMap<bool>* mask_fine, GeographicMap<double>* slope_fine, DInfinityMap* direction_fine, double mindist) {
  if (mask_coarse->getLatitudes() != latitudes || mask_coarse->getLongitudes() != longitudes)
    throw new runtime_error("coarse scale not as expected");

  this->mask_coarse = mask_coarse;
  net = new HydroNet(*mask_coarse);
  out = net->generate(*direction_fine, *mask_fine, *slope_fine, mindist);

  if (stepCallback)
    stepCallback->setup(*this);
}

void SJHydroNetModel::setPrecipitation(PartialConfidenceTemporalGeographicMap<double>* precipitation) {
  this->precipitation = precipitation;
}

void SJHydroNetModel::setTemperature(PartialConfidenceTemporalGeographicMap<double>* surfaceTemp) {
  this->surfaceTemp = surfaceTemp;
}

void SJHydroNetModel::setSnowModel(SnowModel* snowModel) {
  this->snowModel = snowModel;

}

void SJHydroNetModel::setElevation(GeographicMap<double>* elevation) {
  this->elevation = new ScaledGeographicMap<double>(*elevation, latitudes, longitudes, 0.0);
}

void SJHydroNetModel::setPrecipitationScaling(double precipMult) {
  this->precipMult = precipMult;
}

void SJHydroNetModel::setTemperatureAddition(double tempAdd) {
  this->tempAdd = tempAdd;
}

void SJHydroNetModel::setSnowCoverDifference(double snowDiff) {
  this->snowDiff = snowDiff;
}

void SJHydroNetModel::setVerbosity(int verbose) {
  this->verbose = verbose;
}

time_t SJHydroNetModel::getTime() {
  return (time_t) now.getValue();
}

void SJHydroNetModel::runTo(long time) {
  Measure meastime(time, timeind);

  while (now == Measure(0, timeind) || now < meastime) {
    time_t now_time = now.getValue();
    struct tm* ptm = gmtime(&now_time);

    if (now.getValue() != 0) {
      cout << ptm->tm_mday << "/" << ptm->tm_mon+1 << "/" << ptm->tm_year << " < ";

      time_t end_time = meastime.getValue();
      struct tm* end_ptm = gmtime(&end_time);

      cout << end_ptm->tm_mday << "/" << end_ptm->tm_mon+1 << "/" << end_ptm->tm_year << endl;
    }

    stepDay();
  }
}

void SJHydroNetModel::stepDay() {
  if (verbose)
    cout << "Beginning step" << endl;
  if (now.getValue() == 0) {
    // What is the earliest time we can handle?
    cout << "Start time: " << precipitation->getTimes().getMin() << ", " << surfaceTemp->getTimes().getMin() << ", " << snowModel->getTimes().getMin() << endl;
    now = max(max(precipitation->getTimes().getMin(), surfaceTemp->getTimes().getMin()), snowModel->getTimes().getMin());
    cout << "Starting at " << now << endl;
  } else if (timeind == Inds::unixtime)
    now += DividedRange::toTimespan(1);
  else
    now += 1/360.0;

  if (verbose)
    cout << "Rescaling maps" << endl;
  ScaledGeographicMap<double> scaledPrecipitation((*precipitation)[now], latitudes, longitudes, 0.0);
  ScaledGeographicMap<double> scaledSurfaceTemp((*surfaceTemp)[now], latitudes, longitudes, 0.0);
  ScaledGeographicMap<double> scaledSnowCover((*snowModel)[now], latitudes, longitudes, 0.0);

  if (precipMult != 1)
    scaledPrecipitation *= precipMult;
  if (tempAdd != 0)
    scaledSurfaceTemp += tempAdd;

  GeographicMap<double>& fracSnowCover = scaledSnowCover / 100;
  if (snowDiff != 0) {
    fracSnowCover -= snowDiff;
    fracSnowCover *= (fracSnowCover > 0.0);
  }

  /*for (int rr = 0; rr < 13; rr++)
    for (int cc = 0; cc < 32; cc++)
    cout << rr << ", " << cc << ": " << snowModel->debugInfo(rr, cc) << endl;*/

  GeographicMap<double>& fullMeltDegreeDayFactor = (meltDegreeDayFactor + meltDegreeDaySlope * *elevation);

  GeographicMap<double> *newSnowMeltHeightPtr, *newSnowAccumHeightPtr,
    *newRainRunoffHeightPtr, *newMeltRunoffHeightPtr,
    *newRainGroundHeightPtr, *newMeltGroundHeightPtr,
    *newDirectHeightConfPtr;
  stepDayHeight<GeographicMap<double>, GeographicMap<bool> >(scaledPrecipitation, scaledSurfaceTemp,
                                                             fracSnowCover, fullMeltDegreeDayFactor,
                                                             rainRunoffCoefficient, meltRunoffCoefficient,
                                                             groundCoefficient, rainOnSnowCoefficient,
                                                             newSnowMeltHeightPtr, newSnowAccumHeightPtr,
                                                             newRainRunoffHeightPtr, newMeltRunoffHeightPtr,
                                                             newRainGroundHeightPtr, newMeltGroundHeightPtr,
                                                             newDirectHeightConfPtr, verbose);
  GeographicMap<double>& newSnowMeltVolume = *newSnowMeltHeightPtr * *mmdayToVolume;
  GeographicMap<double>& newSnowAccumVolume = *newSnowAccumHeightPtr * *mmdayToVolume;
  GeographicMap<double>& newRainRunoffVolume = *newRainRunoffHeightPtr * *mmdayToVolume;
  GeographicMap<double>& newMeltRunoffVolume = *newMeltRunoffHeightPtr * *mmdayToVolume;
  GeographicMap<double>& newRainGroundVolume = *newRainGroundHeightPtr * *mmdayToVolume;
  GeographicMap<double>& newMeltGroundVolume = *newMeltGroundHeightPtr * *mmdayToVolume;
  GeographicMap<double>& newDirectVolumeConf = *newDirectHeightConfPtr * *mmdayToVolume;

  snowModel->inform(newSnowMeltVolume, newSnowAccumVolume);

  // Extract baseflow
  GeographicMap<double>& newBaseflowRainVolume = *currentGroundRainVolume * groundToBaseflowDay;
  GeographicMap<double>& newBaseflowMeltVolume = *currentGroundMeltVolume * groundToBaseflowDay;

  // Add to ground and extract baseflow
  currentGroundConf = &weightedFraction(newRainGroundVolume + newMeltGroundVolume, newDirectVolumeConf, *currentGroundRainVolume + *currentGroundMeltVolume, *currentGroundConf);
  Transients::abandon(currentGroundConf);

  *currentGroundRainVolume += newRainGroundVolume - newBaseflowRainVolume;
  *currentGroundMeltVolume += newMeltGroundVolume - newBaseflowMeltVolume;

  GeographicMap<double>& newRainVolume = newRainRunoffVolume + newBaseflowRainVolume;
  GeographicMap<double>& newMeltVolume = newMeltRunoffVolume + newBaseflowMeltVolume;
  GeographicMap<double>& newVolumeConf = weightedFraction(newRainRunoffVolume + newMeltRunoffVolume, newDirectVolumeConf, newBaseflowRainVolume + newBaseflowMeltVolume, *currentGroundConf);

  time_t now_time = now.getValue();
  struct tm* ptm = gmtime(&now_time);

  double fracyear = ptm->tm_yday / 365.0;

  if (verbose)
    cout << "Modelling flow" << endl;
  unsigned elapsed = 0;

  net->updateDay(scaledSurfaceTemp, fracSnowCover, newRainVolume, fracyear, surfaceEvaporationFactor, riverEvaporationFactor);

  while (elapsed < DividedRange::toTimespan(1).getValue()) {
    double step = net->calculateMaximumStep();
    if (step < 1)
      step = 1;
    if (step > DividedRange::toTimespan(1).getValue() - elapsed)
      step = DividedRange::toTimespan(1).getValue() - elapsed;

    if (verbose)
      cout << "Elapsed: " << elapsed << ", Timestep: " << step << endl;

    net->step(step, newRainVolume, newMeltVolume, newVolumeConf);

    elapsed += step;
  }

  Transients::clean();

  out->reset();
  outFlowRain.push_back(out->getPrecipVolume());
  outFlowMelt.push_back(out->getMeltVolume());
  if (verbose)
    cout << "Flows: " << out->getPrecipVolume() << ", " << out->getMeltVolume() << endl;

  if (stepCallback)
    stepCallback->post(*this);
}

template <class TTemporal, class TNumeric>
void SJHydroNetModel::runToHeight(TTemporal& precipitation, TTemporal& surfaceTemp,
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
      cout << "Start time: " << precipitation->getTimes().getMin() << ", " << surfaceTemp->getTimes().getMin() << ", " << fracSnowCover->getTimes().getMin() << endl;
      now = max(max(precipitation->getTimes().getMin(), surfaceTemp->getTimes().getMin()), fracSnowCover->getTimes().getMin());
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
    stepDayHeight(nowPrecipitation, nowSurfaceTemp, nowFracSnowCover, fullMeltDegreeDayFactor,
                  rainRunoffCoefficient, meltRunoffCoefficient, groundCoefficient, rainOnSnowCoefficient,
                  newSnowMeltHeightPtr, newSnowAccumHeightPtr, newRainRunoffHeightPtr,
                  newMeltRunoffHeightPtr, newRainGroundHeightPtr, newMeltGroundHeightPtr,
                  newDirectHeightConfPtr);

    newSnowMeltHeight.add(*newSnowMeltHeightPtr);
    newSnowAccumHeight.add(*newSnowAccumHeightPtr);
    newRainRunoffHeight.add(*newRainRunoffHeightPtr);
    newMeltRunoffHeight.add(*newMeltRunoffHeightPtr);
    newRainGroundHeight.add(*newRainGroundHeightPtr);
    newMeltGroundHeight.add(*newMeltGroundHeightPtr);
    newDirectHeightConf.add(*newDirectHeightConfPtr);
  }
}

template <class TNumeric, class TLogical>
void SJHydroNetModel::stepDayHeight(TNumeric& scaledPrecipitation, TNumeric& scaledSurfaceTemp,
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

list<double> SJHydroNetModel::getOutFlowsRain() {
  return outFlowRain;
}

list<double> SJHydroNetModel::getOutFlowsMelt() {
  return outFlowMelt;
}

unsigned SJHydroNetModel::getOutFlowsCount() {
  return outFlowRain.size();
}

list<pair<pair<pair<Measure, Measure>, pair<Measure, Measure> >, pair<bool, double> > > SJHydroNetModel::getAllEdges() {
  list<pair<pair<HydroNode*, HydroNode*>, double> > edges = net->getAllEdges();
  list<pair<pair<pair<Measure, Measure>, pair<Measure, Measure> >, pair<bool, double> > > result;

  list<pair<pair<HydroNode*, HydroNode*>, double> >::iterator it;
  for (it = edges.begin(); it != edges.end(); it++) {
    if (dynamic_cast<HydroOutputNode*>(it->first.second) != NULL && it->first.second != out)
      continue; // ignore node
    pair<Measure, Measure> firstLocation = pair<Measure, Measure>(it->first.first->getLatitude(), it->first.first->getLongitude());
    pair<Measure, Measure> secondLocation = pair<Measure, Measure>(it->first.second->getLatitude(), it->first.second->getLongitude());
    bool isSurface = dynamic_cast<HydroSurfaceNode*>(it->first.first) != NULL;

    result.push_back(pair<pair<pair<Measure, Measure>, pair<Measure, Measure> >, pair<bool, double> >(pair<pair<Measure, Measure>, pair<Measure, Measure> >(firstLocation, secondLocation), pair<bool, double>(isSurface, it->second)));
  }

  return result;
}

// Serializable protocol

/*
istream& operator>>(istream& in, SJHydroNetModel& sink) {
  sink.latitudes = DividedRange::streamExtract(in);
  sink.longitudes = DividedRange::streamExtract(in);
  sink.timeind = Indicator::streamExtract(in);

  in >> *sink.net >> *sink.out;

  // DO NOT INPUT INPUT DATA

  in >> sink.precipMult >> " " >> sink.tempAdd >> " " >> sink.snowDiff;

  // DO NOT INPUT CURRENT STATE
  return in;
}

ostream& operator<<(ostream& os, SJHydroNetModel& source) {
  source.latitudes.streamInsert(os);
  source.longitudes.streamInsert(os);
  source.timeind.streamInsert(os);

  os << *source.net << *source.out;

  // DO NOT OUTPUT INPUT DATA

  os << source.precipMult << " " << source.tempAdd << " " << source.snowDiff;

  // DO NOT OUTPUT CURRENT STATE

  return os;
}
*/

// Helper Functions
GeographicMap<double>& SJHydroNetModel::weightedFraction(GeographicMap<double>& weights1, GeographicMap<double>& values1, GeographicMap<double>& weights2, GeographicMap<double>& values2) {
  return (weights1 * values1 + weights2 * values2).dividedBy(weights1 + weights2, values1);
}
