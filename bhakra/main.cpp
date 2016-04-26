#define USE_DUMMY
#define COMPARE_SNOW

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
#include "BhakraModel.h" // NOTE: haven't redone code for this; decided to lock in USE_DUMMY and COMPARE_SNOW in BhakraModel
#include "../Callbacks.h"

int main(int argc, const char* argv[])
{
  MPI_Init(NULL,NULL); {
    int rank,size;
    MPI_Comm_rank(MCW,&rank);
    MPI_Comm_size(MCW,&size);
    
    SJHydroNetModel* model = makeBhakraModel();

    SJHydroNetModelSaveSomeCells callback("glacierrecord.csv");
    callback.addLocation(31.28, 78.33);
    callback.addLocation(31.4, 78.5);
    callback.addLocation(31.37, 78.49);
    callback.addLocation(30.45, 81.333);
    model->setStepCallback(&callback);

    cout << "Loading bhakra flow" << endl;
    TimeSeries<double>* known = TimeSeries<double>::loadDelimited(DividedRange::withMax(DividedRange::toTime(1963, 1, 1),
                                                                                        DividedRange::toTime(2005, 4, 25),
                                                                                        DividedRange::toTimespan(1).getValue(), Inds::unixtime),
                                                                  "satluj.tsv", NULL, '\t');
    *known *= 2446.56555; // Needs to be multiplied by 2446.57555 for ft^3/s to m^3/day

    cout << "Initialization complete." << endl;

    model->setVerbosity(0);
    
    try {
      model->runTo(DividedRange::toTime(2004, 12, 31));
    } catch (exception& e) {
      cout << "Exception: " << e.what() << endl;
    }
    list<double> predsPrecip = model->getOutFlowsRain();
    list<double> predsMelt = model->getOutFlowsMelt();

    list<double>::iterator it1, it2;
    for (it1 = predsPrecip.begin(), it2 = predsMelt.begin(); it1 != predsPrecip.end() && it2 != predsMelt.end(); it1++, it2++)
      cout << *it1 << ", " << *it2 << endl;
    
  } MPI_Finalize();
}



