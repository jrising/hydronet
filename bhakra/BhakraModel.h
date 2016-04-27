SJHydroNetModel* makeBhakraModel() {
  DividedRange latitudes = DividedRange::withEnds(29.625, 33.875, .25, Inds::lat);
  DividedRange longitudes = DividedRange::withEnds(74.875, 85.125, .25, Inds::lon);

  SJHydroNetModel* model = new SJHydroNetModel(latitudes, longitudes, Inds::unixtime,
					       // meltDegreeDayFactor, meltDegreeDaySlope, rainRunoffCoefficient, meltRunoffCoefficient, groundCoefficient, groundToBaseflowDay, rainOnSnowCoefficient, surfaceEvaporationFactor, riverEvaporationFactor
					       0.633727, 0.000548348, 0.0170707, 3.98606e-06, 0.00510088, 0.0329162, 0.0770918, 0, 0);

  cout << "Setting up model" << endl;
  GeographicMap<float>& slope = *MatrixGeographicMap<float>::loadTIFF(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
								      DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
								      "finalslp.tiff");
  slope /= 1e5; // don't produce transient!

  model->setup(MatrixGeographicMap<double>::loadDelimited(latitudes, longitudes, "mask_new.tsv", NULL, '\t'),
	       MatrixGeographicMap<bool>::loadDelimited(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
							DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
							"mask.tsv", NULL, '\t'),
	       new GeographicMap<double>(slope),
	       new DInfinityMap(*MatrixGeographicMap<float>::loadTIFF(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
								      DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
								      "finalang.tiff")), 10000.0);

  cout << "Loading precipitation" << endl;
  model->setPrecipitation(DelayedPartialConfidenceTemporalGeographicMap<double>::loadDelimited(latitudes, longitudes,
											       DividedRange::withMax(DividedRange::toTime(1948, 1, 1),
														     DividedRange::toTime(2011, 2, 28),
														     DividedRange::toTimespan(1).getValue(), Inds::unixtime),
											       "mergeprecip.tsv",
											       "mergeprecip_conf.tsv", NULL, '\t'));
    
  cout << "Loading temeprature" << endl;
  model->setTemperature(DelayedPartialConfidenceTemporalGeographicMap<double>::loadDelimited(latitudes, longitudes,
											     DividedRange::withMax(DividedRange::toTime(1948, 1, 1),
														   DividedRange::toTime(2011, 2, 8),
														   DividedRange::toTimespan(1).getValue(), Inds::unixtime),
											     "mergetemps.tsv",
											     "mergetemps_conf.tsv", NULL, '\t'));

  cout << "Loading snow cover" << endl;

  DividedRange snowLatitude = DividedRange::withEnds(29.83333, 33.83333, .3333333, Inds::lat);
  DividedRange snowLongitude = DividedRange::withEnds(74.83333, 85.16667, .3333333, Inds::lon);

  ConstantGeographicMap<double> initialVolumes(latitudes, longitudes, 0);

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

  DummyCompareSnowModel* snowCompare = new DummyCompareSnowModel(*snowCover, initialVolumes, snowCoverTime);

  model->setSnowModel(snowCompare);

  cout << "Loading elevation" << endl;
  model->setElevation(MatrixGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
								 DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
								 "elevation.tsv", NULL, '\t'));

  return model;
}
