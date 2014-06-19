#include <iostream>
#include <stdio.h>
#include "../code/SJHydroNetModel.h"
#include "../code/BalanceSnowModel.h"
#include <datastr/GeographicMap.h>
#include <datastr/ConstantGeographicMap.h>
#include <datastr/ScaledGeographicMap.h>
#include <datastr/DelayedTemporalGeographicMap.h>
#include <datastr/DelayedPartialConfidenceTemporalGeographicMap.h>
#include <datastr/TimeSeries.h>
#include <datastr/DividedRange.h>
#include <datastr/FileFormats.h>
#include <tools/hydro/DInfinityMap.h>
#include <metrics/OLS.h>
#include <measure/Inds.h>
#include <measure/Units.h>

using namespace openworld;

int main(int argc, const char* argv[])
{
  MPI_Init(NULL,NULL); {
	int rank,size;
	MPI_Comm_rank(MCW,&rank);
	MPI_Comm_size(MCW,&size);

    Indicator uniformYear = LinearIndicator("Uniform year", Units::s, 0, numeric_limits<double>::max());

    DividedRange precislats = DividedRange::withEnds(29.823100, 34.6906, 0.442500, Inds::lat);
    DividedRange precislons = DividedRange::withEnds(75.797500, 84.6475, 0.442500, Inds::lon);
    DividedRange precistime = DividedRange(1961, 2099, 1.0/360, uniformYear);

    DividedRange latitudes = DividedRange::withEnds(29.625, 33.875, .25, Inds::lat);
    DividedRange longitudes = DividedRange::withEnds(74.875, 85.125, .25, Inds::lon);

    SJHydroNetModel model(latitudes, longitudes, uniformYear,
                          2.98302, -0.000162646, 0.0295818, 0.0117123, 0.00921557, 0.0270544, 0, 0);

    cout << "Setting up model" << endl;
    GeographicMap<float>& slope = *MatrixGeographicMap<float>::loadTIFF(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                                        DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                                        "finalslp.tiff");
    slope /= 1e5; // don't produce transient!

    model.setup(MatrixGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.625, 33.875, .25, Inds::lat),
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
    DelayedTemporalGeographicMap<float>* precip = DelayedTemporalGeographicMap<float>::loadBinary(precislats, precislons, precistime, "/Users/jrising/Dropbox/Himalayan Flooding/PRECIS/precis-q0-himalaya-daily-05216-1961-2098.dat", 4, FileParser<float>::parseNativeEndianFloat);
    model.setPrecipitation(DelayedPartialConfidenceTemporalGeographicMap<double>::convertScale(precip, 86400));
    
    cout << "Loading temeprature" << endl;
    DelayedTemporalGeographicMap<float>* tempmin = DelayedTemporalGeographicMap<float>::loadBinary(precislats, precislons, precistime, "/Users/jrising/Dropbox/Himalayan Flooding/PRECIS/precis-q0-himalaya-daily-03236min-1961-2098.dat", 4, FileParser<float>::parseNativeEndianFloat);
    DelayedTemporalGeographicMap<float>* tempmax = DelayedTemporalGeographicMap<float>::loadBinary(precislats, precislons, precistime, "/Users/jrising/Dropbox/Himalayan Flooding/PRECIS/precis-q0-himalaya-daily-03236max-1961-2098.dat", 4, FileParser<float>::parseNativeEndianFloat);

    DelayedPartialConfidenceTemporalGeographicMap<double>* converted = DelayedPartialConfidenceTemporalGeographicMap<double>::convertSumScale(tempmin, tempmax, 0.5);

    model.setTemperature(converted);

    cout << "Loading snow cover" << endl;
    
    const double initialSnowData[18][42] = {
      // [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42]
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 28890000, 214200000, 365400000, 1.337e+08, 0, 1940000, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18700000, 89200000, 765300000, 43900000, 0, 0, 721400000, 971400000, 848000000, 5.035e+08, 117720000, 0, 0, 0, 0, 77560000},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 500000, 150100000, 144030000, 42530000, 209300000, 453100000, 1001500000, 867600000, 709100000, 1068600000, 814700000, 1075900000, 63470000, 0, 0.000e+00, 0, 0, 12338000, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 480500000, 1325900000, 1152600000, 955700000, 885500000, 744400000, 759200000, 411700000, 0, 0.000e+00, 0, 0, 61200000, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1002900000, 777300000, 3427800000, 763700000, 28300000, 478700000, 166100000, 208190000, 371000000, 177400000, 520600000, 167980000, 602300000, 6.048e+08, 332000000, 304400000, 498800000, 464300000, 557600000, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 34044000, 360730000, 0, 0, 977200000, 9788800000, 134670000, 333600000, 220980000, 0, 0, 0, 0, 0, 31200000, 0, 541800000, 1089800000, 2469400000, 1.480e+09, 386900000, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8720000, 11250000, 579460000, 1194400000, 1528100000, 833800000, 4230200000, 1229500000, 0, 0, 0, 0, 0, 30400000, 770900000, 1288100000, 884100000, 497700000, 271100000, 1021600000, 253200000, 2.459e+08, 28090000, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 20740000, 310000, 121940000, 622000000, 902400000, 675500000, 1606500000, 515370000, 0, 0, 0, 0, 0, 0, 687400000, 807500000, 565700000, 425600000, 189200000, 914400000, 0, 1.387e+08, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100000, 419810000, 955500000, 497200000, 0, 742100000, 399860000, 0, 0, 0, 0, 632000000, 352710000, 485300000, 1102400000, 337400000, 657900000, 532900000, 231700000, 322000000, 0, 0, 0.000e+00, 0, 80820000, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7870000, 694900000, 994500000, 87080000, 36790000, 317950000, 45630000, 0, 0, 0, 0, 14470000, 360800000, 0, 106550000, 0, 323000000, 0, 0, 0, 173380000, 0, 0.000e+00, 0, 0, 0, 0, 120400000, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 190730000, 1312410000, 2181160000, 6500900000, 36920000, 0, 125250000, 9080000, 420020000, 116860000, 1860000, 0, 35890000, 145700000, 0, 332000000, 0, 866500000, 235800000, 289600000, 250400000, 0, 332300000, 2.152e+08, 0, 53100000, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1588760000, 1878509000, 0, 0, 0, 242440000, 18600000, 623700000, 592400000, 671600000, 362600000, 312200000, 476600000, 266500000, 0, 578700000, 253400000, 267000000, 0, 0, 230700000, 0.000e+00, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6526160000, 10474940000, 2910040000, 5425400000, 57490000, 1310000000, 1084900000, 207090000, 0, 0, 0, 0, 0, 0, 0, 485500000, 58500000, 81400000, 73800000, 0, 0, 0, 0.000e+00, 0, 0, 0, 0, 29310000, 20020000},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 64070000, 2875250000, 173930000, 0, 0, 0, 48910000, 0, 0, 0, 1345600000, 86630000, 2223900000, 1925000000, 35900000, 441300000, 270100000, 173900000, 0, 117380000, 0, 0.000e+00, 0, 0, 0, 64410000, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 192240000, 0, 2585120000, 1755054000, 72360000, 0, 0, 0, 116670000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 242370000, 45300000, 0, 0, 414600000, 151930000, 23670000, 0.000e+00, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 485800000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 107990000, 111030000, 0, 0, 371600000, 28424000, 0, 0.000e+00, 0, 0, 0, 0, 0, 0},
      { 0, 80948000, 0, 0, 0, 0, 0, 0, 0, 822130000, 660060000, 905830000, 323030000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 471200000, 503200000, 0, 597900000, 695800000, 588100000, 50800000, 0, 95550000, 0.000e+00, 0, 0, 0, 0, 0, 0}};
    double** initialSnow = new double*[18];
    for (int ii = 0; ii < 18; ii++) {
      initialSnow[ii] = new double[42];
      for (int jj = 0; jj < 42; jj++)
        initialSnow[ii][jj] = initialSnowData[ii][jj];
    }
    MatrixGeographicMap<double> initialVolumes(latitudes, longitudes, initialSnow);

    DividedRange snowLatitude = DividedRange::withEnds(29.83333, 33.83333, .3333333, Inds::lat);
    DividedRange snowLongitude = DividedRange::withEnds(74.83333, 85.16667, .3333333, Inds::lon);
    DelayedTemporalGeographicMap<double>* snowCover = DelayedTemporalGeographicMap<double>::loadDelimited(snowLatitude, snowLongitude,
                                                                                                          DividedRange::withMax(DividedRange::toTime(1988, 1, 1),
                                                       DividedRange::toTime(2005, 5, 1),
                                                       DividedRange::toTimespan(365.25 / 52).getValue(), Inds::unixtime), "snows.tsv", NULL, '\t');

    GeographicMap<double>& firstCover = (*snowCover)[Measure(DividedRange::toTime(1988, 1, 15), Inds::unixtime)];
    ScaledGeographicMap<double> scaledInitialCover(firstCover, latitudes, longitudes, 0.0);

    BalanceSnowModel* snowModel = new BalanceSnowModel(initialVolumes, scaledInitialCover, DividedRange(1988, 2099, 1.0/52.0, uniformYear), 9.178893e-01, 1.396938e-07, -7.527681e-08, 1.119644e-09);
    model.setSnowModel(snowModel);

    cout << "Loading elevation" << endl;
    model.setElevation(MatrixGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                                  DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                                  "elevation.tsv", NULL, '\t'));

    cout << "Initialization complete." << endl;

    try {
      model.runTo(2098);
    } catch (exception& e) {
      cout << "Exception: " << e.what() << endl;
    }
    list<double> predsPrecip = model.getOutFlowsRain(), predsMelt = model.getOutFlowsMelt();

    list<double>::iterator it1, it2;
    for (it1 = predsPrecip.begin(), it2 = predsMelt.begin(); it1 != predsPrecip.end() && it2 != predsMelt.end(); it1++, it2++)
      cout << *it1 << ", " << *it2 << endl;

  } MPI_Finalize();
}
