#include <algorithm>
#include "SJHydroModel.h"
#include <measure/Inds.h>
#include <memory/Transients.h>
#include <utils/Timer.h>

SJHydroModel::SJHydroModel(Indicator timeind, double meltDegreeDayFactor, double meltDegreeDaySlope, double meltDegreeDaySlope, double rainRunoffCoefficient, double meltRunoffCoefficient, double groundCoefficient, double groundToBaseflowDay, double surfaceEvaporationFactor, double riverEvaporationFactor, TNumeric& initialGroundRainVolume, TNumeric& initialGroundMeltVolume, TNumeric& initialGroundConf)
  : timeind(timeind), now(0, timeind), currentGroundRainVolume(initialGroundRainVolume), currentGroundMeltVolume(initialGroundMeltVolume), currentGroundConf(initialGroundConf) {
  precipMult = 1;
  tempAdd = 0;
  snowDiff = 0;

  now = 0;

  this->meltDegreeDayFactor = meltDegreeDayFactor;
  this->meltDegreeDaySlope = meltDegreeDaySlope;
  this->rainRunoffCoefficient = rainRunoffCoefficient;
  this->meltRunoffCoefficient = meltRunoffCoefficient;
  this->groundCoefficient = groundCoefficient;
  this->groundToBaseflowDay = groundToBaseflowDay;
  this->surfaceEvaporationFactor = surfaceEvaporationFactor;
  this->riverEvaporationFactor = riverEvaporationFactor;
}

SJHydroModel::SJHydroModel(SJHydroModel& copy)
  : timeind(copy.timeind), now(copy.now) {
  precipMult = copy.precipMult;
  tempAdd = copy.tempAdd;
  snowDiff = copy.snowDiff;

  meltDegreeDayFactor = copy.meltDegreeDayFactor;
  meltDegreeDaySlope = copy.meltDegreeDaySlope;
  rainRunoffCoefficient = copy.rainRunoffCoefficient;
  meltRunoffCoefficient = copy.meltRunoffCoefficient;
  groundCoefficient = copy.groundCoefficient;
  groundToBaseflowDay = copy.groundToBaseflowDay;
  surfaceEvaporationFactor = copy.surfaceEvaporationFactor;
  riverEvaporationFactor = copy.riverEvaporationFactor;

  outFlowRain = copy.outFlowRain;
  outFlowMelt = copy.outFlowMelt;
}

SJHydroModel::~SJHydroModel() {
}

void SJHydroModel::setPrecipitationScaling(double precipMult) {
  this->precipMult = precipMult;
}

void SJHydroModel::setTemperatureAddition(double tempAdd) {
  this->tempAdd = tempAdd;
}

void SJHydroModel::setSnowCoverDifference(double snowDiff) {
  this->snowDiff = snowDiff;
}

time_t SJHydroModel::getTime() {
  return (time_t) now.getValue();
}

void SJHydroModel::runTo(long time) {
  Measure meastime(time, timeind);

  while (now < meastime) {
    time_t now_time = now.getValue();
    struct tm* ptm = gmtime(&now_time);

    if (now.getValue() != 0)
      cout << ptm->tm_mday << "/" << ptm->tm_mon+1 << "/" << ptm->tm_year << ": " << now << " < " << meastime << endl;

    stepDay();
  }
}

list<double> SJHydroModel::getOutFlowsRain() {
  return outFlowRain;
}

list<double> SJHydroModel::getOutFlowsMelt() {
  return outFlowMelt;
}

unsigned SJHydroModel::getOutFlowsCount() {
  return outFlowRain.size();
}

template<class TNumeric, class TLogical>
SJHydroModel::internalStepDay(TNumeric& nowPrecipitation, TNumeric& nowSurfaceTemp, TNumeric& nowSnowCover) {
  if (precipMult != 1)
    nowPrecipitation *= precipMult;
  if (tempAdd != 0)
    nowSurfaceTemp += tempAdd;

  TNumeric& fracSnowCover = nowSnowCover / 100;
  if (snowDiff != 0) {
    fracSnowCover -= snowDiff;
    fracSnowCover *= (fracSnowCover > 0.0);
  }

  // Add precipitation
  cout << "Determining valid inputs" << endl;
  TLogical& validPrecipitation = nowPrecipitation >= 0.0;
  TLogical& validSurfaceTemp = nowSurfaceTemp >= 100.0;
  TLogical& validSnowCover = nowSnowCover >= 0.0;
  TLogical& aboveZeroCelsius = nowSurfaceTemp >= ZERO_CELSIUS;

  cout << "Adding surface flow" << endl;
  // rain-related calculations
  TNumeric& rainPortion = ((nowSurfaceTemp - SNOW_ALL_TEMPERATURE) / (RAIN_ALL_TEMPERATURE - SNOW_ALL_TEMPERATURE)) * (nowSurfaceTemp > SNOW_ALL_TEMPERATURE) * (nowSurfaceTemp <= RAIN_ALL_TEMPERATURE) + (nowSurfaceTemp > RAIN_ALL_TEMPERATURE);
  TNumeric& newRainAllVolume = rainPortion * nowPrecipitation * *mmdayToVolume * (validPrecipitation * validSurfaceTemp);
  TNumeric& newSnowfreeRainVolume = newRainAllVolume * (1.0 - fracSnowCover) * validSnowCover;
  //TNumeric& newSnowVolume = (1 - rainPortion) * nowPrecipitation * mmdayToVolume * (validPrecipitation * validSurfaceTemp);

  // snow-related calculations
  TNumeric& newRainOnSnowRainVolume = newRainAllVolume * fracSnowCover * validSnowCover; // uses melt runoff, because rain on snow like melt
  TNumeric& newRainOnSnowMeltVolume = (4.2 / 325) * (nowSurfaceTemp - ZERO_CELSIUS) * aboveZeroCelsius * newRainOnSnowRainVolume;
  TNumeric& fullMeltDegreeDayFactor = (meltDegreeDayFactor + meltDegreeDaySlope * *elevation);
  //TNumeric& validMeltDegreeDayFactor = fullMeltDegreeDayFactor * (fullMeltDegreeDayFactor > 0);
  TNumeric& newDegreeDayMeltVolume = fullMeltDegreeDayFactor * fracSnowCover * *mmdayToVolume * (nowSurfaceTemp - ZERO_CELSIUS) * aboveZeroCelsius * validSnowCover;

  TNumeric& newSnowMeltVolume = newDegreeDayMeltVolume + newRainOnSnowMeltVolume;
  TNumeric& newSnowAccumVolume = (1.0 - rainPortion) * nowPrecipitation * *mmdayToVolume * (validPrecipitation * validSurfaceTemp);
  snowModel->inform(newSnowMeltVolume, newSnowAccumVolume);

  // final rain-related and melt-related volumes
  TNumeric& newRainRunoffVolume = rainRunoffCoefficient * newSnowfreeRainVolume + meltRunoffCoefficient * newRainOnSnowRainVolume;
  TNumeric& newMeltRunoffVolume = meltRunoffCoefficient * newSnowMeltVolume;

  TNumeric& newRainGroundVolume = groundCoefficient * ((1 - rainRunoffCoefficient) * newSnowfreeRainVolume + (1 - meltRunoffCoefficient) * newRainOnSnowRainVolume);
  TNumeric& newMeltGroundVolume = groundCoefficient * (1 - meltRunoffCoefficient) * newSnowMeltVolume; // CHANGE from paper: don't use M, which is multiplied by melt coefficient

  TNumeric& newDirectVolumeConf = *(tew_(TNumeric(validPrecipitation * validSurfaceTemp * validSnowCover)));

  // Extract baseflow
  TNumeric& newBaseflowRainVolume = *currentGroundRainVolume * groundToBaseflowDay;
  TNumeric& newBaseflowMeltVolume = *currentGroundMeltVolume * groundToBaseflowDay;

  // Add to ground and extract baseflow
  currentGroundConf = &weightedFraction(newRainGroundVolume + newMeltGroundVolume, newDirectVolumeConf, *currentGroundRainVolume + *currentGroundMeltVolume, *currentGroundConf);
  Transients::abandon(currentGroundConf);

  *currentGroundRainVolume += newRainGroundVolume - newBaseflowRainVolume;
  *currentGroundMeltVolume += newMeltGroundVolume - newBaseflowMeltVolume;

  newRainVolume = newRainRunoffVolume + newBaseflowRainVolume;
  newMeltVolume = newMeltRunoffVolume + newBaseflowMeltVolume;
  newVolumeConf = weightedFraction(newRainRunoffVolume + newMeltRunoffVolume, newDirectVolumeConf, newBaseflowRainVolume + newBaseflowMeltVolume, *currentGroundConf);
}
