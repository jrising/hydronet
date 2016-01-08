//#define USE_DUMMY
//#define COMPARE_SNOW

#include <iostream>
#include <stdio.h>
#include "../SJHydroNetModel.h"
#include "../BackupSnowModel.h"
#ifdef COMPARE_SNOW
#include "../BalanceSnowModel.h"
#include "../CompareSnowModel.h"
#endif
#include <datastr/GeographicMap.h>
#include <datastr/ConstantGeographicMap.h>
#include <datastr/DelayedTemporalGeographicMap.h>
#include <datastr/DelayedPartialConfidenceTemporalGeographicMap.h>
#include <datastr/TimeSeries.h>
#include <datastr/DividedRange.h>
#include <datastr/FileFormats.h>
#include <metrics/OLS.h>
#include <measure/Inds.h>
#ifdef USE_DUMMY
#include "../DummyCompareSnowModel.h"
#endif

int main(int argc, const char* argv[])
{
  MPI_Init(NULL,NULL); {
	int rank,size;
	MPI_Comm_rank(MCW,&rank);
	MPI_Comm_size(MCW,&size);

    DividedRange latitudes = DividedRange::withEnds(29.625, 33.875, .25, Inds::lat);
    DividedRange longitudes = DividedRange::withEnds(74.875, 85.125, .25, Inds::lon);

    SJHydroNetModel model(latitudes, longitudes, Inds::unixtime,
                          // meltDegreeDayFactor, meltDegreeDaySlope, rainRunoffCoefficient, meltRunoffCoefficient, groundCoefficient, groundToBaseflowDay, surfaceEvaporationFactor, riverEvaporationFactor
                          //7.39573, -0.000459404, 0.0263708, 0.00491881, 0.00513267, 0.0478994, 0, 0);
                          3.08558, -0.000108007, 0.028037, 0.0114833, 0.0089908, 0.0281662, 0, 0); // single-year , "allcells.tsv");
                          //2.98302, -0.000162646, 0.0295818, 0.0117123, 0.00921557, 0.0270544, 0, 0); // older single-year

    //model.setSnowCoverDifference(1.0);

    cout << "Setting up model" << endl;
    GeographicMap<float>& slope = *MatrixGeographicMap<float>::loadTIFF(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                                        DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                                        "finalslp.tiff");
    slope /= 1e5; // don't produce transient!

    model.setup(MatrixGeographicMap<double>::loadDelimited(latitudes, longitudes,
                                                           "mask_new.tsv", NULL, '\t'),
                MatrixGeographicMap<bool>::loadDelimited(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                           DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                         "mask.tsv", NULL, '\t'),
                new GeographicMap<double>(slope),
                new DInfinityMap(*MatrixGeographicMap<float>::loadTIFF(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                                       DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                                       "finalang.tiff")), 10000.0);

    /*cout << "Edges:" << endl;
    list<pair<pair<pair<Measure, Measure>, pair<Measure, Measure> >, pair<bool, double> > > edges = model.getAllEdges();
    list<pair<pair<pair<Measure, Measure>, pair<Measure, Measure> >, pair<bool, double> > >::iterator it;
    for (it = edges.begin(); it != edges.end(); it++) {
      cout << "[[" << it->first.first.first.getValue() << ", " << it->first.first.second.getValue() << "], ";
      if (it->first.second.first.getValue() == 0 && it->first.second.second.getValue() == 0)
        cout << "null, ";
      else
        cout << "[" << it->first.second.first.getValue() << ", " << it->first.second.second.getValue() << "], ";
      cout << (it->second.first ? "true, " : "false, ") << it->second.second << "]" << endl;
      }*/

    cout << "Loading precipitation" << endl;
    model.setPrecipitation(DelayedPartialConfidenceTemporalGeographicMap<double>::loadDelimited(latitudes, longitudes,
                                                                                                DividedRange::withMax(DividedRange::toTime(1948, 1, 1),
                                                                                                                      DividedRange::toTime(2011, 2, 28),
                                                                                                                      DividedRange::toTimespan(1).getValue(), Inds::unixtime),
                                                                                                "mergeprecip.tsv",
                                                                                                "mergeprecip_conf.tsv", NULL, '\t'));
    
    cout << "Loading temeprature" << endl;
    model.setTemperature(DelayedPartialConfidenceTemporalGeographicMap<double>::loadDelimited(latitudes, longitudes,
                                                                                              DividedRange::withMax(DividedRange::toTime(1948, 1, 1),
                                                                                                                    DividedRange::toTime(2011, 2, 8),
                                                                                                                    DividedRange::toTimespan(1).getValue(), Inds::unixtime),
                                                                                              "mergetemps.tsv",
                                                                                              "mergetemps_conf.tsv", NULL, '\t'));

    cout << "Loading snow cover" << endl;
    
    DividedRange snowLatitude = DividedRange::withEnds(29.83333, 33.83333, .3333333, Inds::lat);
    DividedRange snowLongitude = DividedRange::withEnds(74.83333, 85.16667, .3333333, Inds::lon);

#ifdef COMPARE_SNOW
#ifndef USE_DUMMY
    const double initialSnowData[18][42] = {
      // [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42]
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 28890000.0, 214200000.0, 365400000.0, 1.337e+08, 0.0, 1940000.0, 0.0, 0.0, 0.0, 0.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18700000.0, 89200000.0, 765300000.0, 43900000.0, 0.0, 0.0, 721400000.0, 971400000.0, 848000000.0, 5.035e+08, 117720000.0, 0.0, 0.0, 0.0, 0.0, 77560000.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 500000.0, 150100000.0, 144030000.0, 42530000.0, 209300000.0, 453100000.0, 1001500000.0, 867600000.0, 709100000.0, 1068600000.0, 814700000.0, 1075900000.0, 63470000.0, 0.0, 0.0, 0.0, 0.0, 12338000.0, 0.0, 0.0, 0.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 480500000.0, 1325900000.0, 1152600000.0, 955700000.0, 885500000.0, 744400000.0, 759200000.0, 411700000.0, 0.0, 0.0, 0.0, 0.0, 61200000.0, 0.0, 0.0, 0.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1002900000.0, 777300000.0, 3427800000.0, 763700000.0, 28300000.0, 478700000.0, 166100000.0, 208190000.0, 371000000.0, 177400000.0, 520600000.0, 167980000.0, 602300000.0, 6.048e+08, 332000000.0, 304400000.0, 498800000.0, 464300000.0, 557600000.0, 0.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 34044000.0, 360730000.0, 0.0, 0.0, 977200000.0, 9788800000.0, 134670000.0, 333600000.0, 220980000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 31200000.0, 0.0, 541800000.0, 1089800000.0, 2469400000.0, 1.480e+09, 386900000.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8720000.0, 11250000.0, 579460000.0, 1194400000.0, 1528100000.0, 833800000.0, 4230200000.0, 1229500000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 30400000.0, 770900000.0, 1288100000.0, 884100000.0, 497700000.0, 271100000.0, 1021600000.0, 253200000.0, 2.459e+08, 28090000.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20740000.0, 310000.0, 121940000.0, 622000000.0, 902400000.0, 675500000.0, 1606500000.0, 515370000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 687400000.0, 807500000.0, 565700000.0, 425600000.0, 189200000.0, 914400000.0, 0.0, 1.387e+08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 419810000.0, 955500000.0, 497200000.0, 0.0, 742100000.0, 399860000.0, 0.0, 0.0, 0.0, 0.0, 632000000.0, 352710000.0, 485300000.0, 1102400000.0, 337400000.0, 657900000.0, 532900000.0, 231700000.0, 322000000.0, 0.0, 0.0, 0.0, 0.0, 80820000.0, 0.0, 0.0, 0.0, 0.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7870000.0, 694900000.0, 994500000.0, 87080000.0, 36790000.0, 317950000.0, 45630000.0, 0.0, 0.0, 0.0, 0.0, 14470000.0, 360800000.0, 0.0, 106550000.0, 0.0, 323000000.0, 0.0, 0.0, 0.0, 173380000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 120400000.0, 0.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 190730000.0, 1312410000.0, 2181160000.0, 6500900000.0, 36920000.0, 0.0, 125250000.0, 9080000.0, 420020000.0, 116860000.0, 1860000.0, 0.0, 35890000.0, 145700000.0, 0.0, 332000000.0, 0.0, 866500000.0, 235800000.0, 289600000.0, 250400000.0, 0.0, 332300000.0, 2.152e+08, 0.0, 53100000.0, 0.0, 0.0, 0.0, 0.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1588760000.0, 1878509000.0, 0.0, 0.0, 0.0, 242440000.0, 18600000.0, 623700000.0, 592400000.0, 671600000.0, 362600000.0, 312200000.0, 476600000.0, 266500000.0, 0.0, 578700000.0, 253400000.0, 267000000.0, 0.0, 0.0, 230700000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6526160000.0, 10474940000.0, 2910040000.0, 5425400000.0, 57490000.0, 1310000000.0, 1084900000.0, 207090000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 485500000.0, 58500000.0, 81400000.0, 73800000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 29310000.0, 20020000.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 64070000.0, 2875250000.0, 173930000.0, 0.0, 0.0, 0.0, 48910000.0, 0.0, 0.0, 0.0, 1345600000.0, 86630000.0, 2223900000.0, 1925000000.0, 35900000.0, 441300000.0, 270100000.0, 173900000.0, 0.0, 117380000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 64410000.0, 0.0, 0.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 192240000.0, 0.0, 2585120000.0, 1755054000.0, 72360000.0, 0.0, 0.0, 0.0, 116670000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 242370000.0, 45300000.0, 0.0, 0.0, 414600000.0, 151930000.0, 23670000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 485800000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 107990000.0, 111030000.0, 0.0, 0.0, 371600000.0, 28424000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      { 0.0, 80948000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 822130000.0, 660060000.0, 905830000.0, 323030000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 471200000.0, 503200000.0, 0.0, 597900000.0, 695800000.0, 588100000.0, 50800000.0, 0.0, 95550000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    double** initialSnow = new double*[18];
    for (int ii = 0; ii < 18; ii++) {
      initialSnow[ii] = new double[42];
      for (int jj = 0; jj < 42; jj++)
        initialSnow[ii][jj] = initialSnowData[ii][jj];
    }
    MatrixGeographicMap<double> initialVolumes(latitudes, longitudes, initialSnow);
#else
    ConstantGeographicMap<double> initialVolumes(latitudes, longitudes, 0);
#endif
#endif

    DividedRange snowCoverTime = DividedRange::withMax(DividedRange::toTime(1988, 1, 1),
                                                       DividedRange::toTime(2003, 5, 1),
                                                       DividedRange::toTimespan(365.25 / 52).getValue(), Inds::unixtime);
    DividedRange fullTime = DividedRange::withMax(DividedRange::toTime(1963, 1, 1),
                                                  DividedRange::toTime(2005, 4, 25),
                                                  DividedRange::toTimespan(1).getValue(), Inds::unixtime);
    BackupSnowModel* snowCover = new BackupSnowModel(DelayedTemporalGeographicMap<double>::loadDelimited(snowLatitude, snowLongitude,
                                                                                  snowCoverTime, "snows.tsv", NULL, '\t'),
                                                     DelayedTemporalGeographicMap<double>::loadDelimited(snowLatitude, snowLongitude,
                                                                                  snowCoverTime, "snows.tsv", NULL, '\t'),
                                                     fullTime,
                                                     .25, Measure(DividedRange::toTime(1988, 1, 15), Inds::unixtime));
#ifndef COMPARE_SNOW
    model.setSnowModel(snowCover);
#else
#ifdef USE_DUMMY
    DummyCompareSnowModel* snowCompare = new DummyCompareSnowModel(*snowCover, initialVolumes, snowCoverTime);
#else
    ScaledGeographicMap<double> scaledInitialCover((*snowCover)[Measure(DividedRange::toTime(1988, 1, 1), Inds::unixtime)], latitudes, longitudes, 0.0);
    BalanceSnowModel* snowModel = new BalanceSnowModel(initialVolumes, scaledInitialCover, snowCoverTime, 9.178893e-01, 1.396938e-07, -7.527681e-08, 1.119644e-09);
    unlink("snowcompare.tsv");
    CompareSnowModel* snowCompare = new CompareSnowModel(*snowCover, *snowModel, snowCoverTime, "snowcompare.tsv");
    model.setSnowModel(snowCompare);
#endif
#endif

    cout << "Loading elevation" << endl;
    model.setElevation(MatrixGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                                  DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                                  "elevation.tsv", NULL, '\t'));

    cout << "Loading bhakra flow" << endl;
    TimeSeries<double>* known = TimeSeries<double>::loadDelimited(DividedRange::withMax(DividedRange::toTime(1963, 1, 1),
                                                                                        DividedRange::toTime(2005, 4, 25),
                                                                                        DividedRange::toTimespan(1).getValue(), Inds::unixtime),
                                                                  "satluj.tsv", NULL, '\t');
    *known *= 2446.56555; // Needs to be multiplied by 2446.57555 for ft^3/s to m^3/day

    cout << "Initialization complete." << endl;

    model.setVerbosity(0);
    
    try {
      model.runTo(DividedRange::toTime(2004, 12, 31));
    } catch (exception& e) {
      cout << "Exception: " << e.what() << endl;
    }
    list<double> predsPrecip = model.getOutFlowsRain(), predsMelt = model.getOutFlowsMelt();

    list<double>::iterator it1, it2;
    for (it1 = predsPrecip.begin(), it2 = predsMelt.begin(); it1 != predsPrecip.end() && it2 != predsMelt.end(); it1++, it2++)
      cout << *it1 << ", " << *it2 << endl;
    
  } MPI_Finalize();
}
