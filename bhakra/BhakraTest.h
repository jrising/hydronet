#ifndef BHAKRA_TEST
#define BHAKRA_TEST

#include <iostream>
#include <stdio.h>
#include <time.h>
#include "../SJHydroNetModel.h"
#include "../BackupSnowModel.h"
#include "../DummyCompareSnowModel.h"
#include <datastr/GeographicMap.h>
#include <datastr/DelayedTemporalGeographicMap.h>
#include <datastr/DelayedPartialConfidenceTemporalGeographicMap.h>
#include <datastr/TimeSeries.h>
#include <datastr/AbstractCollection.h>
#include <datastr/DividedRange.h>
#include <datastr/FileFormats.h>
#include <metrics/OLS.h>
#include <measure/Inds.h>
#include <measure/Measure.h>
#include <dims/Dimensionless.h>
#include <model/ModelTest.h>
#include "BhakraModel.h"

//#define USE_EVAP 1

class BhakraTest : public ModelTest {
 protected:
  SJHydroNetModel* model;
  TimeSeries<double>* known;

  SJHydroNetModel* copy;

 public:
    virtual void calibrate(unsigned maxSteps, string saved = "") {
      calibrateByAsexualEvolutionFullEvaluation(maxSteps, saved);
    }

  virtual map<string, Measure> getParameters() {
    map<string, Measure> params;
    params.insert(pair<string, Measure>("meltDegreeDayFactor", Measure(model->meltDegreeDayFactor, LinearIndicator("degreeDayFactor", Units::none, 0, 10))));
    params.insert(pair<string, Measure>("meltDegreeDaySlope", Measure(model->meltDegreeDaySlope, LinearIndicator("degreeDaySlope", Units::none, -.001, .001))));
    params.insert(pair<string, Measure>("rainRunoffCoefficient", Measure(model->rainRunoffCoefficient, LinearIndicator("runoffCoefficient", Units::none, 0, 1))));
    params.insert(pair<string, Measure>("meltRunoffCoefficient", Measure(model->meltRunoffCoefficient, LinearIndicator("runoffCoefficient", Units::none, 0, 1))));
    params.insert(pair<string, Measure>("groundCoefficient", Measure(model->groundCoefficient, LinearIndicator("groundCoefficient", Units::none, 0, 1))));
    params.insert(pair<string, Measure>("groundToBaseflowDay", Measure(model->groundToBaseflowDay, LinearIndicator("timePortion", Units::none, 0, 1))));
    params.insert(pair<string, Measure>("rainOnSnowCoefficient", Measure(model->rainOnSnowCoefficient, LinearIndicator("runoffCoefficient", Units::none, 0, 1))));

#ifdef USE_EVAP
    params.insert(pair<string, Measure>("surfaceEvaporationFactor", Measure(model->surfaceEvaporationFactor, LinearIndicator("evaporationFactor", Units::none, 0, 2))));
    params.insert(pair<string, Measure>("riverEvaporationFactor", Measure(model->riverEvaporationFactor, LinearIndicator("evaporationFactor", Units::none, 0, 2))));
#endif

    return params;
  }

  virtual void setParameters(map<string, Measure>& params) {
    model->meltDegreeDayFactor = params.find("meltDegreeDayFactor")->second.getValue();
    model->meltDegreeDaySlope = params.find("meltDegreeDaySlope")->second.getValue();
    model->rainRunoffCoefficient = params.find("rainRunoffCoefficient")->second.getValue();
    model->meltRunoffCoefficient = params.find("meltRunoffCoefficient")->second.getValue();
    model->groundCoefficient = params.find("groundCoefficient")->second.getValue();
    model->groundToBaseflowDay = params.find("groundToBaseflowDay")->second.getValue();
    model->rainOnSnowCoefficient = params.find("rainOnSnowCoefficient")->second.getValue();

#ifdef USE_EVAP
    if (params.find("surfaceEvaporationFactor") == params.end()) {
      model->surfaceEvaporationFactor = 0;
      model->riverEvaporationFactor = 0;
    } else {
      model->surfaceEvaporationFactor = params.find("surfaceEvaporationFactor")->second.getValue();
      model->riverEvaporationFactor = params.find("riverEvaporationFactor")->second.getValue();
    }
#endif
  }

  virtual map<string, Measure> suggestParameters(map<string, Measure> params, unsigned count, double distance) {
    map<string, Measure> copy;

#ifdef USE_EVAP
    if (params.find("surfaceEvaporationFactor") == params.end()) {
      copy.insert(pair<string, Measure>("surfaceEvaporationFactor", Measure(.7, LinearIndicator("evaporationFactor", Units::none, 0, 2))));
      copy.insert(pair<string, Measure>("riverEvaporationFactor", Measure(.7, LinearIndicator("evaporationFactor", Units::none, 0, 2))));
    }
#endif

    if (count == 0) {
      copy.insert(pair<string, Measure>("meltDegreeDayFactor", params.find("meltDegreeDayFactor")->second.random()));
      copy.insert(pair<string, Measure>("meltDegreeDaySlope", params.find("meltDegreeDaySlope")->second.random()));
      copy.insert(pair<string, Measure>("rainRunoffCoefficient", params.find("rainRunoffCoefficient")->second.random()));
      copy.insert(pair<string, Measure>("meltRunoffCoefficient", params.find("meltRunoffCoefficient")->second.random()));
      copy.insert(pair<string, Measure>("groundToBaseflowDay", params.find("groundToBaseflowDay")->second.random()));
      copy.insert(pair<string, Measure>("rainOnSnowCoefficient", params.find("rainOnSnowCoefficient")->second.random()));
      copy.insert(pair<string, Measure>("groundCoefficient", params.find("groundCoefficient")->second.random()));
    } else {
      // Modify the parameters
      for (map<string, Measure>::iterator it = params.begin(); it != params.end(); it++) {
        if (rand() % 2 == 1) {
          double mydist = distance;
          if (rand() % 2 == 1) {
            double f = (double) rand() / RAND_MAX;
            mydist = mydist + (1 - mydist) * f;
          }
          Measure attempt = it->second.randomAdjust(mydist);
          cout << "Try " << it->first << " = " << attempt.getValue() << endl;
          copy.insert(pair<string, Measure>(it->first, attempt));
        } else if (rand() % 2 == 1) {
          copy.insert(pair<string, Measure>(it->first, it->second));
        } else {
          copy.insert(pair<string, Measure>(it->first, it->second.random()));
        }
      }
    }

    return copy;
  }

  // Construct HydroNet, evalate to start of Snow data, and save it.
  virtual void prepare() {
    model = makeBhakraModel();

    cout << "Loading bhakra flow" << endl;

    known = TimeSeries<double>::loadDelimited(DividedRange::withMax(DividedRange::toTime(1963, 1, 1),
                                                                    DividedRange::toTime(2005, 4, 25),
                                                                    DividedRange::toTimespan(1).getValue(), Inds::unixtime),
                                              "satluj.tsv", NULL, '\t');
    *known *= 2446.56555; // Needs to be multiplied by 2446.57555 for ft^3/s to m^3/day

    copy = NULL;

    cout << "Initialization complete." << endl;

    model->setVerbosity(0);

    model->runTo(DividedRange::toTime(1988, 1, 1));
  }

  virtual double evaluate() {
    cout << "Copy HydroNet" << endl;
    SJHydroNetModel* copy = new SJHydroNetModel(*model);
    cout << "Using " << copy->meltDegreeDayFactor << ", " << copy->meltDegreeDaySlope << ", " << copy->rainRunoffCoefficient << ", " << copy->meltRunoffCoefficient << ", " << copy->groundCoefficient << ", " << copy->groundToBaseflowDay << ", " << copy->rainOnSnowCoefficient << ", " << copy->surfaceEvaporationFactor << ", " << copy->riverEvaporationFactor << endl;

    copy->setVerbosity(0);

    // Evaluate the glacier points
    double glalats[] = {31.28, 31.4, 31.37, 30.45};
    double glalons[] = {78.33, 78.5, 78.49, 81.333};
    double glasnow[] = {1., 1., 1., 1.};
    MultiMeasure glacierLatitudes(4, glalats, Inds::lat);
    MultiMeasure glacierLongitudes(4, glalons, Inds::lon);
    TemporalAbstractCollection<double>* glacierPrecips = TemporalAbstractCollection<double>::loadMapPoints(copy->getPrecipitation(), glacierLatitudes, glacierLongitudes);
    TemporalAbstractCollection<double>* glacierTemps = TemporalAbstractCollection<double>::loadDelimitedPoints(copy->getTemperature().getTimes(), "glaciertemps.csv");
    TemporalAbstractCollection<double>* glacierCover = TemporalAbstractCollection<double>::loadConstantPoints(copy->getSnowModel().getTimes(), glasnow, 4);
    AbstractCollection<double>* glacierElevation = AbstractCollection<double>::loadMapPoints(copy->getElevation(), glacierLatitudes, glacierLongitudes);

    TemporalAbstractCollection<double> glacierSnowMeltHeight(Inds::unixtime), glacierSnowAccumHeight(Inds::unixtime), glacierRainRunoffHeight(Inds::unixtime),
      glacierMeltRunoffHeight(Inds::unixtime), glacierRainGroundHeight(Inds::unixtime), glacierMeltGroundHeight(Inds::unixtime), glacierDirectHeightConf(Inds::unixtime);

    SJHydroNetModel::runToHeight(*glacierPrecips, *glacierTemps, *glacierCover, *glacierElevation,
                                 copy->meltDegreeDayFactor, copy->meltDegreeDaySlope,
                                 copy->rainRunoffCoefficient, copy->meltRunoffCoefficient,
                                 copy->groundCoefficient, copy->rainOnSnowCoefficient,
                                 glacierSnowMeltHeight,  glacierSnowAccumHeight, glacierRainRunoffHeight,
                                 glacierMeltRunoffHeight, glacierRainGroundHeight, glacierMeltGroundHeight,
                                 glacierDirectHeightConf, Measure(0, Inds::unixtime),
                                 Measure(DividedRange::toTime(2000, 1, 1), Inds::unixtime), true);

    // Accumulate vectors
    TemporalAbstractCollection<double>& snowHeight = (glacierSnowAccumHeight - glacierSnowMeltHeight).cumsum();

    // Compare glacier points
    double p1988 = snowHeight.getSingle(0).get(Measure(DividedRange::toTime(1988, 7, 1), Inds::unixtime));
    double p1989 = snowHeight.getSingle(0).get(Measure(DividedRange::toTime(1989, 7, 1), Inds::unixtime));
    double p1990 = snowHeight.getSingle(0).get(Measure(DividedRange::toTime(1990, 7, 1), Inds::unixtime));
    double p1991 = snowHeight.getSingle(0).get(Measure(DividedRange::toTime(1991, 7, 1), Inds::unixtime));

    double v1988 = -2100;
    double v1989 = -1760;
    double v1990 = -2030;
    double v1991 = -2850;

    cout << "Glacier comparison: " << p1989 - p1988 << ", " << p1990 - p1989 << ", " << p1991 - p1990 <<
      " =?= " << v1989 - v1988 << ", " << v1990 - v1989 << ", " << v1991 - v1990 << endl;

    // Evaluate the gridded model

    unsigned beforeCount = copy->getOutFlowsCount();
    copy->runTo(DividedRange::toTime(2000, 1, 1));

    list<double> predsRain = copy->getOutFlowsRain(), predsMelt = copy->getOutFlowsMelt();
    Matrix<double> preds(copy->getOutFlowsCount() - beforeCount + 3, 1); // 3 more for glacier prediction
    Matrix<double> knownSubset(copy->getOutFlowsCount() - beforeCount + 3, 1); // 3 more for glacier prediction

    list<double>::iterator predsRainIt, predsMeltIt;
    unsigned ii;
    for (predsRainIt = predsRain.begin(), predsMeltIt = predsMelt.begin(), ii = 0;
         predsRainIt != predsRain.end() && predsMeltIt != predsMelt.end(); predsRainIt++, predsMeltIt++, ii++) {
      if (ii < beforeCount)
        continue;
      preds.get(ii - beforeCount, 0) = *predsRainIt + *predsMeltIt;
      knownSubset.get(ii - beforeCount, 0) = known->get(ii);
      cout << known->get(ii) << "\t" << *predsRainIt << "\t" << *predsMeltIt << endl;
    }

    // Add predictions of volume, relative to 1988
    double mm2m3 = 2641759225. / 1000; // m^2 * m / mm

    preds.get(copy->getOutFlowsCount() - beforeCount, 0) = (p1989 - p1988) * mm2m3 / sqrt(365.);
    preds.get(copy->getOutFlowsCount() - beforeCount + 1, 0) = (p1990 - p1989) * mm2m3 / sqrt(365.);
    preds.get(copy->getOutFlowsCount() - beforeCount + 2, 0) = (p1991 - p1990) * mm2m3 / sqrt(365.);

    knownSubset.get(copy->getOutFlowsCount() - beforeCount, 0) = (v1989 - v1988) * mm2m3 / sqrt(365.);
    knownSubset.get(copy->getOutFlowsCount() - beforeCount + 1, 0) = (v1990 - v1989) * mm2m3 / sqrt(365.);
    knownSubset.get(copy->getOutFlowsCount() - beforeCount + 2, 0) = (v1991 - v1990) * mm2m3 / sqrt(365.);

    cout << "Lengths: " << knownSubset.getRows() << " vs. " << preds.getRows() << endl;
    double rsqr = OLS::calcRSqr(knownSubset, preds);
    cout << "RSqr: " << rsqr << endl;

    delete copy;

    return rsqr;
  }

  virtual double partialEvaluate() {
    if (!copy)
      copy = new SJHydroNetModel(*model);
    else {
      copy->meltDegreeDayFactor = model->meltDegreeDayFactor;
      copy->meltDegreeDaySlope = model->meltDegreeDaySlope;
      copy->rainRunoffCoefficient = model->rainRunoffCoefficient;
      copy->meltRunoffCoefficient = model->meltRunoffCoefficient;
      copy->groundCoefficient = model->groundCoefficient;
      copy->groundToBaseflowDay = model->groundToBaseflowDay;
      copy->rainOnSnowCoefficient = model->rainOnSnowCoefficient;

#ifdef USE_EVAP
      copy->surfaceEvaporationFactor = model->surfaceEvaporationFactor;
      copy->riverEvaporationFactor = model->riverEvaporationFactor;
#endif
    }

    cout << "Using " << copy->meltDegreeDayFactor << ", " << copy->meltDegreeDaySlope << ", " << copy->rainRunoffCoefficient << ", " << copy->meltRunoffCoefficient << ", " << copy->groundCoefficient << ", " << copy->groundToBaseflowDay << ", " << copy->rainOnSnowCoefficient << ", " << copy->surfaceEvaporationFactor << ", " << copy->riverEvaporationFactor << endl;

    double rsqr2 = yearEvaluate();

    time_t now = copy->getTime();
    struct tm* tm = gmtime(&now);
    if (tm->tm_year >= 100) {
      delete copy;
      copy = NULL;
    }

    return rsqr2;
  }

  virtual double yearEvaluate() {
    unsigned beforeCount = copy->getOutFlowsCount();

    time_t now = copy->getTime();
    struct tm* tm = gmtime(&now);
    if (tm->tm_year == 90) {
      // drop 1990 - 1991, with incomplete cover data
      model->runTo(DividedRange::toTime(1992, 1, 1));
      beforeCount = copy->getOutFlowsCount();
      tm->tm_year = 93;
    } else
      tm->tm_year++;

    copy->runTo(mktime(tm));

    list<double> predsRain = copy->getOutFlowsRain(), predsMelt = copy->getOutFlowsMelt();
    Matrix<double> preds(copy->getOutFlowsCount() - beforeCount, 1);

    list<double>::iterator predsRainIt, predsMeltIt;
    unsigned ii;
    for (predsRainIt = predsRain.begin(), predsMeltIt = predsMelt.begin(), ii = 0;
         predsRainIt != predsRain.end() && predsMeltIt != predsMelt.end(); predsRainIt++, predsMeltIt++, ii++) {
      if (ii < beforeCount)
        continue;
      preds.get(ii - beforeCount, 0) = *predsRainIt + *predsMeltIt;
      cout << known->get(ii) << "\t" << *predsRainIt << "\t" << *predsMeltIt << endl;
    }

    Matrix<double>& knownSubset = known->subset(beforeCount, 0, preds.getRows(), 1);
    double rsqr = OLS::calcRSqr(knownSubset, preds);
    cout << "RSqr: " << rsqr << endl;

    return rsqr;
  }
};

#endif
