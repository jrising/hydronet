#include <algorithm>
#include "SJHydroNetModel.h"
#include <measure/Inds.h>
#include <memory/Transients.h>
#include <utils/Timer.h>
#include <datastr/ConstantGeographicMap.h>

SJHydroNetModel::SJHydroNetModel(DividedRange latitudes, DividedRange longitudes, Indicator timeind, double meltDegreeDayFactor, double meltDegreeDaySlope, double rainRunoffCoefficient, double meltRunoffCoefficient, double groundCoefficient, double groundToBaseflowDay, double surfaceEvaporationFactor, double riverEvaporationFactor, string allcells)
  : latitudes(latitudes), longitudes(longitudes), timeind(timeind), now(0, timeind) {
  precipitation = NULL;
  surfaceTemp = NULL;
  snowModel = NULL;
  elevation = NULL;
  mask_coarse = NULL;

  precipMult = 1;
  tempAdd = 0;
  snowDiff = 0;

  now = 0;
  this->allcells = allcells;

  this->meltDegreeDayFactor = meltDegreeDayFactor;
  this->meltDegreeDaySlope = meltDegreeDaySlope;
  this->rainRunoffCoefficient = rainRunoffCoefficient;
  this->meltRunoffCoefficient = meltRunoffCoefficient;
  this->groundCoefficient = groundCoefficient;
  this->groundToBaseflowDay = groundToBaseflowDay;
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
  out = (HydroOutputNode*) HydroNode::getCopy(copy.out, translate);

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
  allcells = copy.allcells;

  meltDegreeDayFactor = copy.meltDegreeDayFactor;
  meltDegreeDaySlope = copy.meltDegreeDaySlope;
  rainRunoffCoefficient = copy.rainRunoffCoefficient;
  meltRunoffCoefficient = copy.meltRunoffCoefficient;
  groundCoefficient = copy.groundCoefficient;
  groundToBaseflowDay = copy.groundToBaseflowDay;
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

  if (!allcells.empty()) {
    ofstream cellfile;
    cellfile.open(allcells.c_str(), ios::out);
    for (unsigned rr = 1; rr < latitudes.count() - 1; rr++)
      for (unsigned cc = 1; cc < longitudes.count() - 1; cc++)
        cellfile << latitudes.getCellCenter(rr).getValue() << "\t" << longitudes.getCellCenter(cc).getValue() << "\t";
    cellfile << endl;
    cellfile.close();
  }
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

time_t SJHydroNetModel::getTime() {
  return (time_t) now.getValue();
}

void SJHydroNetModel::runTo(long time) {
  Measure meastime(time, timeind);

  while (now < meastime) {
    time_t now_time = now.getValue();
    struct tm* ptm = gmtime(&now_time);

    if (now.getValue() != 0)
      cout << ptm->tm_mday << "/" << ptm->tm_mon+1 << "/" << ptm->tm_year << ": " << now << " < " << meastime << endl;

    stepDay();
  }
}

void SJHydroNetModel::stepDay() {
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

  // Add precipitation
  cout << "Determining valid inputs" << endl;
  GeographicMap<bool>& validPrecipitation = scaledPrecipitation >= 0.0;
  GeographicMap<bool>& validSurfaceTemp = scaledSurfaceTemp >= 100.0;
  GeographicMap<bool>& validSnowCover = scaledSnowCover >= 0.0;
  GeographicMap<bool>& aboveZeroCelsius = scaledSurfaceTemp >= ZERO_CELSIUS;

  cout << "Adding surface flow" << endl;
  // rain-related calculations
  GeographicMap<double>& rainPortion = ((scaledSurfaceTemp - SNOW_ALL_TEMPERATURE) / (RAIN_ALL_TEMPERATURE - SNOW_ALL_TEMPERATURE)) * (scaledSurfaceTemp > SNOW_ALL_TEMPERATURE) * (scaledSurfaceTemp <= RAIN_ALL_TEMPERATURE) + (scaledSurfaceTemp > RAIN_ALL_TEMPERATURE);
  GeographicMap<double>& newRainAllVolume = rainPortion * scaledPrecipitation * *mmdayToVolume * (validPrecipitation * validSurfaceTemp);
  GeographicMap<double>& newSnowfreeRainVolume = newRainAllVolume * (1.0 - fracSnowCover) * validSnowCover;
  //GeographicMap<double>& newSnowVolume = (1 - rainPortion) * scaledPrecipitation * mmdayToVolume * (validPrecipitation * validSurfaceTemp);

  // snow-related calculations
  GeographicMap<double>& newRainOnSnowRainVolume = newRainAllVolume * fracSnowCover * validSnowCover; // uses melt runoff, because rain on snow like melt
  GeographicMap<double>& newRainOnSnowMeltVolume = (4.2 / 325) * (scaledSurfaceTemp - ZERO_CELSIUS) * aboveZeroCelsius * newRainOnSnowRainVolume;
  GeographicMap<double>& fullMeltDegreeDayFactor = (meltDegreeDayFactor + meltDegreeDaySlope * *elevation);
  //GeographicMap<double>& validMeltDegreeDayFactor = fullMeltDegreeDayFactor * (fullMeltDegreeDayFactor > 0);
  GeographicMap<double>& newDegreeDayMeltVolume = fullMeltDegreeDayFactor * fracSnowCover * *mmdayToVolume * (scaledSurfaceTemp - ZERO_CELSIUS) * aboveZeroCelsius * validSnowCover;

  //cout << "A: 1, 29: " << scaledPrecipitation.getCellConst(1, 29) << ", " << scaledSurfaceTemp.getCellConst(1, 29) << ", " << newDegreeDayMeltVolume.getCellConst(1, 29) << ", " << newRainOnSnowMeltVolume.getCellConst(1, 29) << ", 1-" << rainPortion.getCellConst(1, 29) << " * " << scaledPrecipitation.getCellConst(1, 29) << " * " << mmdayToVolume->getCellConst(1, 29) << " * " << (validPrecipitation.getCellConst(1, 29) * validSurfaceTemp.getCellConst(1, 29)) << " = " << (1.0 - rainPortion.getCellConst(1, 29)) * scaledPrecipitation.getCellConst(1, 29) * mmdayToVolume->getCellConst(1, 29) * (validPrecipitation.getCellConst(1, 29) * validSurfaceTemp.getCellConst(1, 29)) << endl;

  GeographicMap<double>& newSnowMeltVolume = newDegreeDayMeltVolume + newRainOnSnowMeltVolume;
  GeographicMap<double>& newSnowAccumVolume = (1.0 - rainPortion) * scaledPrecipitation * *mmdayToVolume * (validPrecipitation * validSurfaceTemp);
  snowModel->inform(newSnowMeltVolume, newSnowAccumVolume);

  // final rain-related and melt-related volumes
  GeographicMap<double>& newRainRunoffVolume = rainRunoffCoefficient * newSnowfreeRainVolume + meltRunoffCoefficient * newRainOnSnowRainVolume;
  GeographicMap<double>& newMeltRunoffVolume = meltRunoffCoefficient * newSnowMeltVolume;

  GeographicMap<double>& newRainGroundVolume = groundCoefficient * ((1 - rainRunoffCoefficient) * newSnowfreeRainVolume + (1 - meltRunoffCoefficient) * newRainOnSnowRainVolume);
  GeographicMap<double>& newMeltGroundVolume = groundCoefficient * (1 - meltRunoffCoefficient) * newSnowMeltVolume; // CHANGE from paper: don't use M, which is multiplied by melt coefficient

  GeographicMap<double>& newDirectVolumeConf = *(tew_(GeographicMap<double>(validPrecipitation * validSurfaceTemp * validSnowCover)));

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

  cout << "Modelling flow" << endl;
  unsigned elapsed = 0;

  net->updateDay(scaledSurfaceTemp, fracSnowCover, newRainVolume, fracyear, surfaceEvaporationFactor, riverEvaporationFactor);

  while (elapsed < DividedRange::toTimespan(1).getValue()) {
    double step = net->calculateMaximumStep();
    if (step < 1)
      step = 1;
    if (step > DividedRange::toTimespan(1).getValue() - elapsed)
      step = DividedRange::toTimespan(1).getValue() - elapsed;

    cout << "Elapsed: " << elapsed << ", Timestep: " << step << endl;

    net->step(step, newRainVolume, newMeltVolume, newVolumeConf);

    elapsed += step;
  }

  Transients::clean();

  out->reset();
  outFlowRain.push_back(out->getPrecipVolume());
  outFlowMelt.push_back(out->getMeltVolume());
  cout << "Flows: " << out->getPrecipVolume() << ", " << out->getMeltVolume() << endl;
  
  if (!allcells.empty()) {
    ofstream cellfile;
    cellfile.open(allcells.c_str(), ios::out | ios::app);
    for (unsigned rr = 1; rr < latitudes.count() - 1; rr++)
      for (unsigned cc = 1; cc < longitudes.count() - 1; cc++) {
        double precipVolume, meltVolume;
        net->sumNodeMapVolumes(rr, cc, precipVolume, meltVolume);
      
        cellfile << precipVolume << "\t" << meltVolume << "\t";
      }
    cellfile << endl;
    cellfile.close();
  }
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

// Helper Functions
GeographicMap<double>& SJHydroNetModel::weightedFraction(GeographicMap<double>& weights1, GeographicMap<double>& values1, GeographicMap<double>& weights2, GeographicMap<double>& values2) {
  return (weights1 * values1 + weights2 * values2).dividedBy(weights1 + weights2, values1);
}
