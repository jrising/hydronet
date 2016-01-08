#ifndef BHAKRA_TEST
#define BHAKRA_TEST

#include <iostream>
#include <stdio.h>
#include <time.h>
#include "../SJHydroNetModel.h"
#include "../BackupSnowModel.h"
#include <datastr/GeographicMap.h>
#include <datastr/DelayedTemporalGeographicMap.h>
#include <datastr/DelayedPartialConfidenceTemporalGeographicMap.h>
#include <datastr/TimeSeries.h>
#include <datastr/DividedRange.h>
#include <datastr/FileFormats.h>
#include <metrics/OLS.h>
#include <measure/Inds.h>
#include <measure/Measure.h>
#include <dims/Dimensionless.h>
#include <model/ModelTest.h>

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

    /*if (count == 0) {
      copy.insert(pair<string, Measure>("meltDegreeDayFactor", params.find("meltDegreeDayFactor")->second.random()));
      copy.insert(pair<string, Measure>("meltDegreeDaySlope", params.find("meltDegreeDaySlope")->second.random()));
      copy.insert(pair<string, Measure>("rainRunoffCoefficient", params.find("rainRunoffCoefficient")->second.random()));
      copy.insert(pair<string, Measure>("meltRunoffCoefficient", params.find("meltRunoffCoefficient")->second.random()));
      copy.insert(pair<string, Measure>("groundToBaseflowDay", params.find("groundToBaseflowDay")->second.random()));
      copy.insert(pair<string, Measure>("rainOnSnowCoefficient", params.find("rainOnSnowCoefficient")->second.random()));
      copy.insert(pair<string, Measure>("groundCoefficient", params.find("groundCoefficient")->second.random()));
      } else {*/
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
        } else if (rand() % 2 == 1)
          copy.insert(pair<string, Measure>(it->first, it->second));
        else
          copy.insert(pair<string, Measure>(it->first, it->second.random()));	  
      }
      //}

    return copy;
  }

  // Construct HydroNet, evalate to start of Snow data, and save it.
  virtual void prepare() {
    model = new SJHydroNetModel(DividedRange::withEnds(29.625, 33.875, .25, Inds::lat),
                                DividedRange::withEnds(74.875, 85.125, .25, Inds::lon), Inds::unixtime,
                                3.06469, 1.0025e-05, 0.0305725, 0.0159555, 0.0104312, 0.0136191, (4.2 / 325), 0, 0);
                                //2.81594, 0.0, 0.0331202, 0.0166743, 0.0105677, 0.0153526, 0.0, 0.0);

    cout << "Setting up model" << endl;
    GeographicMap<float>& slope = *MatrixGeographicMap<float>::loadTIFF(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                                        DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                                        "finalslp.tiff");
    slope /= 1e5; // don't produce transient!

    model->setup(MatrixGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.625, 33.875, .25, Inds::lat),
                                                           DividedRange::withEnds(74.875, 85.125, .25, Inds::lon),
                                                           "mask_new.tsv", NULL, '\t'),
                MatrixGeographicMap<bool>::loadDelimited(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                           DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                         "mask.tsv", NULL, '\t'),
                new GeographicMap<double>(slope),
                new DInfinityMap(*MatrixGeographicMap<float>::loadTIFF(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                                       DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                                       "finalang.tiff")), 10000.0);

    cout << "Loading precipitation" << endl;
    model->setPrecipitation(DelayedPartialConfidenceTemporalGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.625, 33.875, .25, Inds::lat),
                                                                                                DividedRange::withEnds(74.875, 85.125, .25, Inds::lon),
                                                                                                DividedRange::withMax(DividedRange::toTime(1948, 1, 1),
                                                                                                                      DividedRange::toTime(2011, 2, 28),
                                                                                                                      DividedRange::toTimespan(1).getValue(), Inds::unixtime),
                                                                                                "mergeprecip.tsv",
                                                                                                "mergeprecip_conf.tsv", NULL, '\t'));
    
    cout << "Loading temeprature" << endl;
    model->setTemperature(DelayedPartialConfidenceTemporalGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.625, 33.875, .25, Inds::lat),
                                                                                              DividedRange::withEnds(74.875, 85.125, .25, Inds::lon),
                                                                                              DividedRange::withMax(DividedRange::toTime(1948, 1, 1),
                                                                                                                    DividedRange::toTime(2011, 2, 8),
                                                                                                                    DividedRange::toTimespan(1).getValue(), Inds::unixtime),
                                                                                              "mergetemps.tsv",
                                                                                              "mergetemps_conf.tsv", NULL, '\t'));

    cout << "Loading snow cover" << endl;
    BackupSnowModel* snowCover = new BackupSnowModel(DelayedTemporalGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.83333, 33.83333, .3333333, Inds::lat),
                                                                                                         DividedRange::withEnds(74.83333, 85.16667, .3333333, Inds::lon),
                                                                                                         DividedRange::withMax(DividedRange::toTime(1988, 1, 1),
                                                                                                                               DividedRange::toTime(2003, 5, 1),
                                                                                                                               DividedRange::toTimespan(365.25 / 52).getValue(), Inds::unixtime),
                                                                                                         "snows.tsv", NULL, '\t'),
                                                     DelayedTemporalGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.83333, 33.83333, .3333333, Inds::lat),
                                                                                                         DividedRange::withEnds(74.83333, 85.16667, .3333333, Inds::lon),
                                                                                                         DividedRange::withMax(DividedRange::toTime(1988, 1, 1),
                                                                                                                               DividedRange::toTime(2003, 5, 1),
                                                                                                                               DividedRange::toTimespan(365.25 / 52).getValue(), Inds::unixtime),
                                                                                                         "snows.tsv", NULL, '\t'),
                                                     DividedRange::withMax(DividedRange::toTime(1963, 1, 1),
                                                                           DividedRange::toTime(2005, 4, 25),
                                                                           DividedRange::toTimespan(1).getValue(), Inds::unixtime),
                                                     .25, Measure(DividedRange::toTime(1988, 1, 15), Inds::unixtime));
    model->setSnowModel(snowCover);

    cout << "Loading elevation" << endl;
    model->setElevation(MatrixGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                                   DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                                   "elevation.tsv", NULL, '\t'));

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

    unsigned beforeCount = copy->getOutFlowsCount();
    copy->runTo(DividedRange::toTime(2000, 1, 1));

    list<double> predsRain = copy->getOutFlowsRain(), predsMelt = copy->getOutFlowsMelt();
    Matrix<double> preds(copy->getOutFlowsCount() - beforeCount, 1);

    list<double>::iterator predsRainIt, predsMeltIt;
    unsigned ii;
    for (predsRainIt = predsRain.begin(), predsMeltIt = predsMelt.begin(), ii = 0;
         predsRainIt != predsRain.end() && predsMeltIt != predsMelt.end(); predsRainIt++, predsMeltIt++, ii++) {
      if (ii < beforeCount)
	continue;
      preds.get(ii - beforeCount, 0) = *predsRainIt + *predsMeltIt;
      cout << known->get(ii, 0) << "\t" << *predsRainIt << "\t" << *predsMeltIt << endl;
    }

    Matrix<double>& knownSubset = known->subset(beforeCount, 0, preds.getRows(), 1);
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
      cout << known->get(ii, 0) << "\t" << *predsRainIt << "\t" << *predsMeltIt << endl;
    }
    
    Matrix<double>& knownSubset = known->subset(beforeCount, 0, preds.getRows(), 1);
    double rsqr = OLS::calcRSqr(knownSubset, preds);
    cout << "RSqr: " << rsqr << endl;

    return rsqr;
  }
};

#endif
