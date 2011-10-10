#include <algorithm>
#include "HydroModel.h"
#include <indicator/Inds.h>
#include <memory/Transients.h>

HydroModel::HydroModel(DividedRange latitudes, DividedRange longitudes, double precipCoefficient, double meltCoefficient, double latFlow, double lonFlow)
  : latitudes(latitudes), longitudes(longitudes), now(0, Inds::unixtime) {
  riverVolume = new MatrixGeographicMap<double>(latitudes, longitudes);
  surfaceVolume = new MatrixGeographicMap<double>(latitudes, longitudes);
  riverConf = new MatrixGeographicMap<double>(latitudes, longitudes);
  surfaceConf = new MatrixGeographicMap<double>(latitudes, longitudes);

  riverVolume->loadConstantInto(0);
  surfaceVolume->loadConstantInto(0);
  riverConf->loadConstantInto(0);
  surfaceConf->loadConstantInto(0);
    
  elevation = NULL;
  slope = NULL;
  direction = NULL;

  precipitation = NULL;
  surfaceTemp = NULL;
  snowCover = NULL;

  this->precipCoefficient = precipCoefficient;
  this->meltCoefficient = meltCoefficient;

  now = 0;

  this->rrFlow = latitudes.inRange(latFlow);
  this->ccFlow = longitudes.inRange(lonFlow);

  minSurface = .1;
  minRiver = .001;
  maxSurfaceVelocity = 1;
  maxRiverVelocity = 5;
}

HydroModel::HydroModel(DividedRange latitudes, DividedRange longitudes, double precipCoefficient, double meltCoefficient, unsigned rrFlow, unsigned ccFlow)
  : latitudes(latitudes), longitudes(longitudes), now(0, Inds::unixtime) {
  riverVolume = new MatrixGeographicMap<double>(latitudes, longitudes);
  surfaceVolume = new MatrixGeographicMap<double>(latitudes, longitudes);
  riverConf = new MatrixGeographicMap<double>(latitudes, longitudes);
  surfaceConf = new MatrixGeographicMap<double>(latitudes, longitudes);

  riverVolume->loadConstantInto(0);
  surfaceVolume->loadConstantInto(0);
  riverConf->loadConstantInto(0);
  surfaceConf->loadConstantInto(0);
    
  elevation = NULL;
  slope = NULL;
  direction = NULL;

  precipitation = NULL;
  surfaceTemp = NULL;
  snowCover = NULL;

  this->precipCoefficient = precipCoefficient;
  this->meltCoefficient = meltCoefficient;

  now = 0;

  this->rrFlow = rrFlow;
  this->ccFlow = ccFlow;

  minSurface = .1;
  minRiver = .001;
  maxSurfaceVelocity = 1;
  maxRiverVelocity = 5;
}

void HydroModel::setElevation(GeographicMap<double>* elevation) {
  this->elevation = elevation;
}

void HydroModel::setDInfinity(GeographicMap<double>* slope, GeographicMap<double>* direction) {
  this->direction = direction;
  cout << "Direction from " << direction->min() << " to " << direction->max() << endl;
  this->slope = slope;
}

void HydroModel::setPrecipitation(PartialConfidenceTemporalGeographicMap<double>* precipitation) {
  this->precipitation = precipitation;
}

void HydroModel::setTemperature(PartialConfidenceTemporalGeographicMap<double>* surfaceTemp) {
  this->surfaceTemp = surfaceTemp;
}

void HydroModel::setSnowCover(TemporalGeographicMap<double>* snowCover) {
  this->snowCover = snowCover;
}

void HydroModel::runTo(time_t time) {
  Measure unixtime(time, Inds::unixtime);

  while (now < unixtime) {
    cout << now << " < " << unixtime << endl;
    stepDay();
  }
}

void HydroModel::stepDay() {
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

  double mindist = calcMinimumDistance(*riverVolume);
  cout << "Minimum distance: " << mindist << endl;

  //cout << "Slopes:" << endl << *slope << endl;
  //cout << "Direction:" << endl << *direction << endl;

  cout << "Modelling flow" << endl;
  unsigned elapsed = 0;

  double outFlowDay = 0;
  while (elapsed < DividedRange::toTimespan(1)) {
    //cout << "Elapsed: " << elapsed << ", Flew: " << outFlowDay << endl;

    // determine how far we can go
    GeographicMap<double>& riverVelocity = calcManningRiver(riverVolume); // in m/s
    //cout << "River Velocity:" << endl << riverVelocity << endl;

    GeographicMap<double>& surfaceVelocity = calcManningSurface(surfaceVolume); // in m/s
    //cout << "Surface Velocity:" << endl << surfaceVelocity << endl;

    //cout << "Volume: " << riverVolume->getCellConst(rrFlow, ccFlow) << " + " << surfaceVolume->getCellConst(rrFlow, ccFlow) << " at " << 
    //riverVelocity.getCellConst(rrFlow, ccFlow) << ", " << surfaceVelocity.getCellConst(rrFlow, ccFlow) << endl;

    //cout << "Max velocites: " << riverVelocity.max() << ", " << surfaceVelocity.max() << endl;
    double step = mindist / max(riverVelocity.max(), surfaceVelocity.max()); // in s
    if (step > DividedRange::toTimespan(1) - elapsed)
      step = DividedRange::toTimespan(1) - elapsed;
    if (step < 1)
      step = 1;

    // move the liquid
    MatrixGeographicMap<double>* riverVolumeAfter = new MatrixGeographicMap<double>(latitudes, longitudes);
    riverVolumeAfter->loadConstantInto(0);
    MatrixGeographicMap<double>* surfaceVolumeAfter = new MatrixGeographicMap<double>(latitudes, longitudes);
    MatrixGeographicMap<double>* riverConfAfter = new MatrixGeographicMap<double>(latitudes, longitudes);
    MatrixGeographicMap<double>* riverConfWeightAfter = new MatrixGeographicMap<double>(latitudes, longitudes);
    riverConfAfter->loadConstantInto(0);
    riverConfWeightAfter->loadConstantInto(0);

    for (unsigned rr = 0; rr < latitudes.count(); rr++)
      for (unsigned cc = 0; cc < longitudes.count(); cc++) {
        if (rr == 0 || cc == 0 || rr == latitudes.count() - 1 || cc == longitudes.count() - 1) {
          riverVolumeAfter->getCell(rr, cc) = 0;
          surfaceVolumeAfter->getCell(rr, cc) = 0;
          continue;
        }

        double surface = surfaceVolume->getCellConst(rr, cc), river = riverVolume->getCellConst(rr, cc);
        if (surface < minSurface && river < minRiver) {
          riverVolumeAfter->getCell(rr, cc) += river;
          riverConfAfter->getCell(rr, cc) += riverConf->getCellConst(rr, cc);
          riverConfWeightAfter->getCell(rr, cc) += 1;
          surfaceVolumeAfter->getCell(rr, cc) = surface;
        } else {
          unsigned rr0, cc0, rr1, cc1; // rr0, cc0 is always straight, rr1, cc1 is diagonal
          double portion0; // portion1 = 1 - portion0

          double dir = direction->getCellConst(rr, cc);
          if (dir < 0 || dir > 2 * PI)
            dir = 2 * PI * ((double) rand() / (double) RAND_MAX);
          
          if (dir < M_PI / 4) {
            rr0 = rr;
            cc0 = cc1 = cc + 1;
            rr1 = rr - 1;
            portion0 = 1 - dir / (M_PI / 4);
          } else if (dir < M_PI / 2) {
            rr0 = rr1 = rr - 1;
            cc0 = cc;
            cc1 = cc + 1;
            portion0 = (dir - M_PI / 4) / (M_PI / 4);
          } else if (dir < 3 * M_PI / 4) {
            rr0 = rr1 = rr - 1;
            cc0 = cc;
            cc1 = cc - 1;
            portion0 = 1 - (dir - M_PI / 2) / (M_PI / 4);
          } else if (dir < M_PI) {
            rr0 = rr;
            cc0 = cc1 = cc - 1;
            rr1 = rr - 1;
            portion0 = (dir - 3 * M_PI / 4) / (M_PI / 4);
          } else if (dir < 5 * M_PI / 4) {
            rr0 = rr;
            cc0 = cc1 = cc - 1;
            rr1 = rr + 1;
            portion0 = 1 - (dir - M_PI) / (M_PI / 4);
          } else if (dir < 3 * M_PI / 2) {
            rr0 = rr1 = rr + 1;
            cc0 = cc;
            cc1 = cc - 1;
            portion0 = (dir - 5 * M_PI / 4) / (M_PI / 4);
          } else if (dir < 7 * M_PI / 4) {
            rr0 = rr1 = rr + 1;
            cc0 = cc;
            cc1 = cc + 1;
            portion0 = 1 - (dir - 3 * M_PI / 2) / (M_PI / 4);
          } else {
            rr0 = rr;
            cc0 = cc1 = cc + 1;
            rr1 = rr + 1;
            portion0 = (dir - 7 * M_PI / 4) / (M_PI / 4);
          }
        
          double distanceStraight = riverVolume->calcDistance(rr, cc, rr0, cc0);
          double distanceDiagonal = riverVolume->calcDistance(rr, cc, rr1, cc1);
          
          if (river >= minRiver) {
            double rivervel = riverVelocity.getCellConst(rr, cc);
            double riverPortionStraight = (rivervel * step) / distanceStraight;
            double riverPortionDiagonal = (rivervel * step) / distanceDiagonal;
            riverVolumeAfter->getCell(rr0, cc0) += riverPortionStraight * portion0 * river;
            riverVolumeAfter->getCell(rr1, cc1) += riverPortionDiagonal * (1 - portion0) * river;
            riverVolumeAfter->getCell(rr, cc) += (1 - riverPortionStraight * portion0 - riverPortionDiagonal * (1 - portion0)) * river;
            riverConfAfter->getCell(rr0, cc0) += riverPortionStraight * portion0 * riverConf->getCellConst(rr, cc);
            riverConfAfter->getCell(rr1, cc1) += riverPortionDiagonal * (1 - portion0) * riverConf->getCellConst(rr, cc);
            riverConfAfter->getCell(rr, cc) += (1 - riverPortionStraight * portion0 - riverPortionDiagonal * (1 - portion0)) * riverConf->getCellConst(rr, cc);
            riverConfWeightAfter->getCell(rr0, cc0) += riverPortionStraight * portion0;
            riverConfWeightAfter->getCell(rr1, cc1) += riverPortionDiagonal * (1 - portion0);
            riverConfWeightAfter->getCell(rr, cc) += (1 - riverPortionStraight * portion0 - riverPortionDiagonal * (1 - portion0));

            if (rr == rrFlow && cc == ccFlow)
              outFlowDay += (riverPortionStraight * portion0 + riverPortionDiagonal * (1 - portion0)) * river;
          } else {
            riverVolumeAfter->getCell(rr, cc) += river;
            riverConfAfter->getCell(rr, cc) += riverConf->getCellConst(rr, cc);
            riverConfWeightAfter->getCell(rr, cc) += 1;
          }

          if (surface >= minSurface) {
            double surfacevel = surfaceVelocity.getCellConst(rr, cc);
            double surfacePortionStraight = (surfacevel * step) / distanceStraight;
            double surfacePortionDiagonal = (surfacevel * step) / distanceDiagonal;
            riverVolumeAfter->getCell(rr0, cc0) += surfacePortionStraight * portion0 * surface;
            riverVolumeAfter->getCell(rr1, cc1) += surfacePortionDiagonal * (1 - portion0) * surface;
            surfaceVolumeAfter->getCell(rr, cc) = (1 - surfacePortionStraight * portion0 - surfacePortionDiagonal * (1 - portion0)) * surface;
            riverConfAfter->getCell(rr0, cc0) += surfacePortionStraight * portion0 * surfaceConf->getCellConst(rr, cc);
            riverConfAfter->getCell(rr1, cc1) += surfacePortionDiagonal * (1 - portion0) * surfaceConf->getCellConst(rr, cc);
            riverConfWeightAfter->getCell(rr0, cc0) += surfacePortionStraight * portion0;
            riverConfWeightAfter->getCell(rr1, cc1) += surfacePortionDiagonal * (1 - portion0);
          
            if (rr == rrFlow && cc == ccFlow)
              outFlowDay += (surfacePortionStraight * portion0 + surfacePortionDiagonal * (1 - portion0)) * surface;
          } else
            surfaceVolumeAfter->getCell(rr, cc) = surface;
        }

        double stepFraction = step / (double) DividedRange::toTimespan(1);
        double addition = newVolume.getCellConst(rr, cc) * stepFraction;
        if (addition > 0) {
          double riverArea = calcRiverWidth(riverVolume->getCellConst(rr, cc)) * 1000;
          double cellArea = riverVolume->calcArea(rr, cc);
          double riverPortion = min(riverArea / cellArea, 1.0);

          riverVolumeAfter->getCell(rr, cc) += addition * riverPortion;
          riverConfAfter->getCell(rr, cc) += newVolumeConf.getCellConst(rr, cc) * stepFraction;
          riverConfWeightAfter->getCell(rr, cc) += stepFraction;

          if (riverPortion < 1) {
            surfaceVolumeAfter->getCell(rr, cc) += addition * (1 - riverPortion);
            surfaceConf->getCell(rr, cc) = (surfaceConf->getCellConst(rr, cc) * (2 - stepFraction) + newVolumeConf.getCellConst(rr, cc) * stepFraction) / 2;
          }
        }
      }

    *riverConfAfter /= *riverConfWeightAfter;

    delete riverVolume;
    delete surfaceVolume;
    delete riverConf;
    riverVolume = riverVolumeAfter;
    surfaceVolume = surfaceVolumeAfter;
    riverConf = riverConfAfter;
    delete riverConfWeightAfter;

    elapsed += step;

    cout << "Cleaning memory, elapsed: " << elapsed << endl;
    Transients::clean();
  }

  outFlow.push_back(outFlowDay);

  //cout << "River Volume:" << endl << *riverVolume << endl;
  //cout << "Surface Volume:" << endl << *surfaceVolume << endl;
}

list<double> HydroModel::getOutFlows() {
  return outFlow;
}

GeographicMap<double>* HydroModel::getRiverVolume() {
  return riverVolume;
}

// Utility functions

// Slope is drop / run
double HydroModel::calcManning(double coeff, double radius, double slope) {
  if (radius <= 0 || slope <= 0)
    return 0.0;
  double vel = (1 / coeff) * pow(radius, 2.0/3.0) * sqrt(slope);
  //if (vel > maxRiverVelocity)
  //cout << "1/c (" << radius << ")^2/3 (" << slope << ")^1/2 = " << vel << endl;

  return vel;
}

double HydroModel::calcRiverWidth(double volume) {
  return 2 * sqrt((3 / M_PI) * volume / 1000);
}

GeographicMap<double>& HydroModel::calcManningRiver(GeographicMap<double>* volume) {
  GeographicMap<double>* result = tew_(MatrixGeographicMap<double>(volume->getLatitudes(), volume->getLongitudes()));

  for (unsigned rr = 0; rr < volume->getLatitudes().count(); rr++)
    for (unsigned cc = 0; cc < volume->getLongitudes().count(); cc++) {
      if (rr == 0 || cc == 0 || rr == latitudes.count() - 1 || cc == longitudes.count() - 1) {
        result->getCell(rr, cc) = 0;
        continue;
      }

      double vol = volume->getCellConst(rr, cc);
      if (vol < minRiver) {
        result->getCell(rr, cc) = 0;
        continue;
      }

      // 1/3 of circle calculation
      double width = calcRiverWidth(vol);
      if (width <= 0) {
        result->getCell(rr, cc) = 0;
        continue;
      }

      double r, height;
      if (width > 1000) {
        r = vol / volume->calcArea(rr, cc);
        height = r;
      } else {
        r = (vol / 1000) / (M_PI * width / 3);
        height = width / 2;
      }

      double vel = calcManning(.033, r, max(height / (2 * 1000), slope->getCellConst(rr, cc)));
      result->getCell(rr, cc) = min(vel, maxRiverVelocity); //90 * sqrt(height)); // terminal velocity
    }

  return *result;
}

GeographicMap<double>& HydroModel::calcManningSurface(GeographicMap<double>* volume) {
  GeographicMap<double>* result = tew_(MatrixGeographicMap<double>(volume->getLatitudes(), volume->getLongitudes()));

  for (unsigned rr = 0; rr < volume->getLatitudes().count(); rr++)
    for (unsigned cc = 0; cc < volume->getLongitudes().count(); cc++) {
      if (rr == 0 || cc == 0 || rr == latitudes.count() - 1 || cc == longitudes.count() - 1) {
        result->getCell(rr, cc) = 0;
        continue;
      }

      double vol = volume->getCellConst(rr, cc);
      if (vol < minSurface) {
        result->getCell(rr, cc) = 0;
        continue;
      }
      
      double height = vol / volume->calcArea(rr, cc);
      double vel = calcManning(.025, height, max(height / (2 * 1000), slope->getCellConst(rr, cc)));
      result->getCell(rr, cc) = min(vel, maxSurfaceVelocity); //90 * sqrt(height / 2)); //(water terminal velocity)
    }

  return *result;
}

double HydroModel::calcMinimumDistance(GeographicMap<double>& map) {
  double mindist = map.calcDistance(0, 0, 1, 0);

  for (unsigned rr = 1; rr < map.getLatitudes().count() - 1; rr++)
    mindist = min(mindist, map.calcDistance(rr, 0, rr + 1, 0));

  for (unsigned cc = 0; cc < map.getLongitudes().count() - 1; cc++)
    mindist = min(mindist, map.calcDistance(0, cc, 0, cc + 1));  

  return mindist;
}
