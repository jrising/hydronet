#include <iostream>
#include <stdio.h>
#include <datastr/TimeSeries.h>
#include <datastr/DividedRange.h>
#include <datastr/FileFormats.h>
#include <measure/Inds.h>
#include <datastr/GeographicMap.h>
#include <datastr/DelayedPartialConfidenceTemporalGeographicMap.h>
#include <datastr/PartialConfidenceTemporalGeographicMap.h>

using namespace openworld;

int main(int argc, const char* argv[])
{
  const char* precipsOut = argv[1];
  const char* tempsOut = argv[2];

  double lat0 = atof(argv[3]);
  double lon0 = atof(argv[4]);

  DividedRange latitudes = DividedRange::withEnds(29.625, 33.875, .25, Inds::lat);
  DividedRange longitudes = DividedRange::withEnds(74.875, 85.125, .25, Inds::lon);

  cout << "Loading precipitation" << endl;
  PartialConfidenceTemporalGeographicMap<double>* precipitation =
    DelayedPartialConfidenceTemporalGeographicMap<double>::loadDelimited(latitudes, longitudes,
                                                                                                DividedRange::withMax(DividedRange::toTime(1901, 1, 1),
                                                                                                                      DividedRange::toTime(2010, 12, 31), // 2011, 1, 1
                                                                                                                      DividedRange::toTimespan(1).getValue(), Inds::unixtime),
                                                                           "../bhakra/mergeprecip.tsv",
                                                                           "../bhakra/mergeprecip_conf.tsv", NULL, '\t');

  cout << "Loading temperature" << endl;
  PartialConfidenceTemporalGeographicMap<double>* surfaceTemp = DelayedPartialConfidenceTemporalGeographicMap<double>::loadDelimited(latitudes, longitudes,
                                                                                              DividedRange::withMax(DividedRange::toTime(1948, 1, 1),
                                                                                                                    DividedRange::toTime(2011, 2, 8),
                                                                                                                    DividedRange::toTimespan(1).getValue(), Inds::unixtime),
                                                                                              "../bhakra/mergetemps.tsv",
                                                                                              "../bhakra/mergetemps_conf.tsv", NULL, '\t');

  Measure lat0m(lat0, Inds::lat);
  Measure lon0m(lon0, Inds::lon);

  cout << "Extracting... " << lat0m << ", " << lon0m << endl;
  TimeSeries<double>* precips = precipitation->getTimeSeries(lat0m, lon0m);
  TimeSeries<double>* temps = surfaceTemp->getTimeSeries(lat0m, lon0m);
  cout << "Saving..." << endl;
  precips->saveDelimited(precipsOut);
  temps->saveDelimited(tempsOut);
}

