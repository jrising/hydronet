#include <iostream>
#include <stdio.h>
#include "../SJHydroNetModel.h"
#include "../BackupSnowModel.h"
#include "../BalanceSnowModel.h"
#include "../CompareSnowModel.h"
#include <datastr/TimeSeries.h>
#include <datastr/DividedRange.h>
#include <datastr/FileFormats.h>
#include <metrics/OLS.h>
#include <measure/Inds.h>

int main(int argc, const char* argv[])
{
  SJHydroPointModel model(Inds::unixtime,
                          // meltDegreeDayFactor, meltDegreeDaySlope, rainRunoffCoefficient, meltRunoffCoefficient, groundCoefficient, groundToBaseflowDay, surfaceEvaporationFactor, riverEvaporationFactor
                          7.39573, -0.000459404, 0.0263708, 0.00491881, 0.00513267, 0.0478994, 0, 0);

  cout << "Loading precipitation" << endl;
  model.setPrecipitation(TimeSeries<double>::loadDelimited(DividedRange::withMax(DividedRange::toTime(1901, 1, 1),
                                                                                 DividedRange::toTime(2010, 1, 1),
                                                                                 DividedRange::toTimespan(1).getValue(), Inds::unixtime),
                                                           "mergeprecip.tsv",
                                                           "mergeprecip_conf.tsv", NULL, '\t'));

  cout << "Loading temeprature" << endl;
  model.setTemperature(TimeSeries<double>::loadDelimited(DividedRange::withMax(DividedRange::toTime(1948, 1, 1),
                                                                               DividedRange::toTime(2011, 2, 8),
                                                                               DividedRange::toTimespan(1).getValue(), Inds::unixtime),
                                                         "mergetemps.tsv",
                                                         "mergetemps_conf.tsv", NULL, '\t'));

  cout << "Loading snow cover" << endl;

  DividedRange snowCoverTime = DividedRange::withMax(DividedRange::toTime(1988, 1, 1),
                                                     DividedRange::toTime(2003, 5, 1),
                                                     DividedRange::toTimespan(365.25 / 52).getValue(), Inds::unixtime);
  DividedRange fullTime = DividedRange::withMax(DividedRange::toTime(1963, 1, 1),
                                                DividedRange::toTime(2005, 4, 25),
                                                DividedRange::toTimespan(1).getValue(), Inds::unixtime);
  BackupSnowModel* snowCover = new BackupSnowModel(TimeSeries<double>::loadDelimited(snowLatitude, snowLongitude,
                                                                                  snowCoverTime, "snows.tsv", NULL, '\t'),
                                                     DelayedTimeSeries<double>::loadDelimited(snowLatitude, snowLongitude,
                                                                                  snowCoverTime, "snows.tsv", NULL, '\t'),
                                                     fullTime,
                                                     .25, Measure(DividedRange::toTime(1988, 1, 15), Inds::unixtime));
    model.setSnowModel(snowCover);

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
