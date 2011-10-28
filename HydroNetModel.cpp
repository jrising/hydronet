#include <algorithm>
#include "HydroNetModel.h"
#include <indicator/Inds.h>
#include <memory/Transients.h>
#include <utils/Timer.h>

HydroNetModel::HydroNetModel(DividedRange latitudes, DividedRange longitudes, double precipCoefficient, double meltCoefficient)
  : latitudes(latitudes), longitudes(longitudes), now(0, Inds::unixtime) {    
  precipitation = NULL;
  surfaceTemp = NULL;
  snowCover = NULL;

  this->precipCoefficient = precipCoefficient;
  this->meltCoefficient = meltCoefficient;

  now = 0;
}

void HydroNetModel::setup(GeographicMap<double>* mask_coarse, GeographicMap<bool>* mask_fine, GeographicMap<double>* slope_fine, DInfinityMap* direction_fine) {
  if (mask_coarse.getLatitudes() != latitudes || mask_coarse.getLongitudes() != longitudes)
    throw new runtime_error("coarse scale not as expected");

  net = new HydroNet(*mask);
  out = net->generate(direction_fine, mask_fine, slope_fine)
}

void HydroNetModel::setPrecipitation(PartialConfidenceTemporalGeographicMap<double>* precipitation) {
  this->precipitation = precipitation;
}

void HydroNetModel::setTemperature(PartialConfidenceTemporalGeographicMap<double>* surfaceTemp) {
  this->surfaceTemp = surfaceTemp;
}

void HydroNetModel::setSnowCover(TemporalGeographicMap<double>* snowCover) {
  this->snowCover = snowCover;
}

void HydroNetModel::runTo(time_t time) {
  Measure unixtime(time, Inds::unixtime);

  while (now < unixtime) {
    cout << now << " < " << unixtime << endl;
    stepDay();
  }
}

void HydroNetModel::stepDay() {
  cout << "Beginning step" << endl;
  if (now.getValue() == 0) {
    now = max(max(precipitation->getTimes().getMin(), surfaceTemp->getTimes().getMin()), snowCover->getTimes().getMin());
    cout << "Starting at " << now << endl;
  } else
    now += DividedRange::toTimespan(1);

  cout << "Rescaling maps" << endl;
  ScaledGeographicMap<double> scaledPrecipitation((*precipitation)[now], latitudes, longitudes);
  ScaledGeographicMap<double> scaledSurfaceTemp((*surfaceTemp)[now], latitudes, longitudes);
  ScaledGeographicMap<double> scaledSnowCover((*snowCover)[now], latitudes, longitudes);

  // Add precipitation
  cout << "Determining valid inputs" << endl;
  GeographicMap<bool>& validPrecipitation = scaledPrecipitation >= 0.0;
  GeographicMap<bool>& validSurfaceTemp = scaledSurfaceTemp >= 100.0;
  GeographicMap<bool>& validSnowCover = scaledSnowCover >= 0.0;

  cout << "Adding surface flow" << endl;
  MatrixGeographicMap<double> newVolume(latitudes, longitudes);
  MatrixGeographicMap<double> newVolumeConf(latitudes, longitudes);

  for (unsigned rr = 1; rr < latitudes.count() - 1; rr++)
    for (unsigned cc = 1; cc < longitudes.count() - 1; cc++) {
      double precip = precipCoefficient * scaledPrecipitation.getCellConst(rr, cc) * (scaledSurfaceTemp.getCellConst(rr, cc) >= 272.15) * (validPrecipitation.getCellConst(rr, cc) * validSurfaceTemp.getCellConst(rr, cc));
      double melt = meltCoefficient * scaledSnowCover.getCellConst(rr, cc) * (scaledSurfaceTemp.getCellConst(rr, cc) - 273.15) * (validSurfaceTemp.getCellConst(rr, cc) * validSnowCover.getCellConst(rr, cc));
      newVolume.getCell(rr, cc) = precip + melt;

      double conf = (validPrecipitation.getCellConst(rr, cc) * validSurfaceTemp.getCellConst(rr, cc) + validSurfaceTemp.getCellConst(rr, cc) * validSnowCover.getCellConst(rr, cc)) / 2;
      newVolumeConf.getCell(rr, cc) = conf;
    }

  cout << "Modelling flow" << endl;
  unsigned elapsed = 0;

  double outFlowDay = 0;
  while (elapsed < DividedRange::toTimespan(1)) {
    double step = calculateMaximumStep();
    if (step < 1)
      step = 1;
    if (step > DividedRange::toTimespan(1) - elapsed)
      step = DividedRange::toTimespan(1) - elapsed;

    cout << "Elapsed: " << elapsed << ", Timestep: " << step << endl;

    net->step(dt, step / (double) DividedRange::toTimespan(1), newVolume, newVolumeConf);

    elapsed += step;

    Transients::clean();
  }

  outFlow.push_back(outFlowDay);
}

list<double> HydroNetModel::getOutFlows() {
  return outFlow;
}
