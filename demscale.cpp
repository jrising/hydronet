#include "DInfinityMap.h"

using namespace openworld;

// call as demscale elevation.tsv elevation_new.tsv
int main(int argc, const char* argv[]) {
  MPI_Init(NULL,NULL); {
	int rank,size;
	MPI_Comm_rank(MCW,&rank);
	MPI_Comm_size(MCW,&size);

    GeographicMap<double>& elevation =
      *MatrixGeographicMap<double>::loadDelimited(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                  DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                  argv[1], NULL, '\t');
    MatrixGeographicMap<short> elevation_new(DividedRange::withEnds(29.625, 33.875, .25, Inds::lat),
                                             DividedRange::withEnds(74.875, 85.125, .25, Inds::lon));

    for (unsigned rr_new = 0; rr_new < elevation_new.getLatitudes().count(); rr_new++)
      for (unsigned cc_new = 0; cc_new < elevation_new.getLongitudes().count(); cc_new++) {
        double total = 0, count = 0;
        // loop through elements within this frame
        for (Measure lat = elevation_new.getLatitudes().getCellMin(rr_new) + elevation.getLatitudes().getWidths() / 2;
             lat < elevation_new.getLatitudes().getCellMax(rr_new); lat += elevation.getLatitudes().getWidths())
          for (Measure lon = elevation_new.getLongitudes().getCellMin(cc_new) + elevation.getLongitudes().getWidths() / 2;
               lon < elevation_new.getLongitudes().getCellMax(cc_new); lon += elevation.getLongitudes().getWidths()) {
            if (lat < elevation.getLatitudes().getMin() || lat >= elevation.getLatitudes().getMax() ||
                lon < elevation.getLongitudes().getMin() || lon >= elevation.getLongitudes().getMax())
              continue;

            total += elevation.getDouble(lat, lon);
            count++;
          }
        
        elevation_new.getCell(rr_new, cc_new) = total / count;
      }

    elevation_new.getValues().saveDelimited(argv[2], FileFormatter<short>::formatSimple, '\t');
  } MPI_Finalize();
}

void checkBidirectional(unsigned rr0, unsigned cc0, unsigned rr1, unsigned cc1, DInfinityMap& direction_new) {
  double amount = direction_new.flowsIntoAmount(rr0, cc0, rr1, cc1) * direction_new.flowsIntoAmount(rr1, cc1, rr0, cc0);
  if (amount > .1)
    cout << "Bidirectional at " << rr0 << ", " << cc0 << " and " << rr1 << ", " << cc1 << ": " << amount << endl;
}
