#include <algorithm>
#include "SJHydroPointModel.h"
#include <measure/Inds.h>
#include <memory/Transients.h>
#include <utils/Timer.h>

SJHydroPointModel::SJHydroPointModel(Indicator timeind, double meltDegreeDayFactor, double meltDegreeDaySlope, double rainRunoffCoefficient, double meltRunoffCoefficient, double groundCoefficient, double groundToBaseflowDay, double surfaceEvaporationFactor, double riverEvaporationFactor, double mmdayToVolume)
  : SJHydroModel(timeind, meltDegreeDayFactor, meltDegreeDaySlope, rainRunoffCoefficient, meltRunoffCoefficient, groundCoefficient, groundToBaseflowDay, surfaceEvaporationFactor, riverEvaporationFactor) {
  precipitation = NULL;
  surfaceTemp = NULL;
  snowModel = NULL;
  elevation = 0;

  currentGroundRainVolume = 0;
  currentGroundMeltVolume = 0;
  currentGroundConf = 0;

  this->mmdayToVolume = mmdayToVolume;
}

SJHydroPointModel::SJHydroPointModel(SJHydroPointModel& copy)
  : SJHydroModel(copy) {

  precipitation = copy.precipitation->clone();
  surfaceTemp = copy.surfaceTemp->clone();
  snowModel = copy.snowModel->clone();
  elevation = copy.elevation;

  mmdayToVolume = copy.mmdayToVolume;

  currentGroundRainVolume = copy.currentGroundRainVolume;
  currentGroundMeltVolume = copy.currentGroundMeltVolume;
  currentGroundConf = copy.currentGroundConf;

  outFlowRain = copy.outFlowRain;
  outFlowMelt = copy.outFlowMelt;
}

SJHydroPointModel::~SJHydroPointModel() {
  delete precipitation;
  delete surfaceTemp;
  delete snowModel;
}

void SJHydroPointModel::setPrecipitation(PartialConfidenceTimeSeries<double>* precipitation) {
  this->precipitation = precipitation;
}

void SJHydroPointModel::setTemperature(PartialConfidenceTimeSeries<double>* surfaceTemp) {
  this->surfaceTemp = surfaceTemp;
}

void SJHydroPointModel::setSnowModel(SnowPointModel* snowModel) {
  this->snowModel = snowModel;
}

void SJHydroPointModel::stepDay() {
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
  double nowPrecipitation = precipitation[now];
  double nowSurfaceTemp = surfaceTemp[now];
  double nowSnowCover = snowModel[now];

  internalStepDay(nowPrecipitation, nowSurfaceTemp, nowSnowCover);

  time_t now_time = now.getValue();
  struct tm* ptm = gmtime(&now_time);
  
  double fracyear = ptm->tm_yday / 365.0;

  outFlowRain.push_back(nowRainVolume);
  outFlowMelt.push_back(nowMeltVolume);
  cout << "Flows: " << nowRainVolume << ", " << nowMeltVolume << endl;

  Transients::clean();
}

list<double> SJHydroPointModel::getOutFlowsRain() {
  return outFlowRain;
}

list<double> SJHydroPointModel::getOutFlowsMelt() {
  return outFlowMelt;
}

unsigned SJHydroPointModel::getOutFlowsCount() {
  return outFlowRain.size();
}
