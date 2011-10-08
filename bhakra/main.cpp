#include <iostream>
#include <stdio.h>
#include "HydroModel.h"
#include <datastr/GeographicMap.h>
#include <datastr/DelayedTemporalGeographicMap.h>
#include <datastr/DelayedPartialConfidenceTemporalGeographicMap.h>
#include <datastr/TimeSeries.h>
#include <datastr/DividedRange.h>
#include <datastr/FileFormats.h>
#include <metrics/OLS.h>
#include <indicator/Inds.h>

int main(int argc, const char* argv[])
{
  MPI_Init(NULL,NULL); {
	int rank,size;
	MPI_Comm_rank(MCW,&rank);
	MPI_Comm_size(MCW,&size);

    // XXX: Did I check that all my sources of precip and temp are scaled the same?
    // XXX: Check that I get a reasonable value when I scale precip as it should be

    HydroModel model(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                     DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                     0.05111236083067, 0.04331253217195, 31.399803, 76.367394);
    cout << "Loading elevation" << endl;
    model.setElevation(MatrixGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                                  DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                                  "../regress/elevation.tsv", NULL, '\t'));

    cout << "Loading d-infinity" << endl;
    GeographicMap<float>& slope = *MatrixGeographicMap<float>::loadTIFF(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                                        DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                                        "../regress2/finalslp.tiff");
    slope /= 1e5; // don't produce transient!

    model.setDInfinity(new GeographicMap<double>(slope),
                       new GeographicMap<double>(*MatrixGeographicMap<float>::loadTIFF(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                                                      DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                                                      "../regress2/finalang.tiff")));

    cout << "Loading precipitation" << endl;
    model.setPrecipitation(DelayedPartialConfidenceTemporalGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.625, 33.875, .25, Inds::lat),
                                                                                                DividedRange::withEnds(74.875, 85.125, .25, Inds::lon),
                                                                                                DividedRange::withMax(DividedRange::toTime(1901, 1, 1),
                                                                                                                      DividedRange::toTime(2011, 2, 28),
                                                                                                                      DividedRange::toTimespan(1), Inds::unixtime),
                                                                                                "../regress2/mergeprecip.tsv",
                                                                                                "../regress2/mergeprecip_conf.tsv", NULL, '\t'));
    
    cout << "Loading temeprature" << endl;
    model.setTemperature(DelayedPartialConfidenceTemporalGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.625, 33.875, .25, Inds::lat),
                                                                                              DividedRange::withEnds(74.875, 85.125, .25, Inds::lon),
                                                                                              DividedRange::withMax(DividedRange::toTime(1948, 1, 1),
                                                                                                                    DividedRange::toTime(2011, 2, 8),
                                                                                                                    DividedRange::toTimespan(1), Inds::unixtime),
                                                                                              "../regress2/mergetemps.tsv",
                                                                                              "../regress2/mergetemps_conf.tsv", NULL, '\t'));

    cout << "Loading snow cover" << endl;
    model.setSnowCover(DelayedTemporalGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.83333, 33.83333, .3333333, Inds::lat),
                                                                           DividedRange::withEnds(74.83333, 85.16667, .3333333, Inds::lon),
                                                                           DividedRange::withMax(DividedRange::toTime(1988, 1, 1),
                                                                                                 DividedRange::toTime(2003, 5, 20),
                                                                                                 DividedRange::toTimespan(7), Inds::unixtime),
                                                                           "../regress/snows.tsv", NULL, '\t'));

    cout << "Loading bhakra flow" << endl;
    TimeSeries<double>* known = TimeSeries<double>::loadDelimited(DividedRange::withMax(DividedRange::toTime(1963, 1, 1),
                                                                                        DividedRange::toTime(2005, 4, 25),
                                                                                        DividedRange::toTimespan(1), Inds::unixtime),
                                                                  "../regress2/bhakra.tsv", NULL, '\t');
    cout << "Initialization complete." << endl;

    model.runTo(DividedRange::toTime(1988, 2, 1));
    list<double> preds = model.getOutFlows();

    for (list<double>::iterator it = preds.begin(); it != preds.end(); it++)
      cout << *it << endl;

    //printf("RSqr: %f\n", OLS::calcRSqr(known, preds));
    
  } MPI_Finalize();
}
