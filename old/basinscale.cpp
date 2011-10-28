#include "DInfinityMap.h"

using namespace openworld;

void checkBidirectional(unsigned rr0, unsigned cc0, unsigned rr1, unsigned cc1, DInfinityMap& direction_new);

// call as basinscale ang.tiff mask.tsv ang_new.tsv mask_new.tsv
int main(int argc, const char* argv[]) {
  MPI_Init(NULL,NULL); {
	int rank,size;
	MPI_Comm_rank(MCW,&rank);
	MPI_Comm_size(MCW,&size);

    DInfinityMap direction(*MatrixGeographicMap<float>::loadTIFF(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                                 DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                                 argv[1]), false);
    GeographicMap<bool>& mask = *MatrixGeographicMap<bool>::loadDelimited(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                                          DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                                          argv[2], NULL, '\t');

    MatrixGeographicMap<double> direction_new(DividedRange::withEnds(29.625, 33.875, .25, Inds::lat),
                                              DividedRange::withEnds(74.875, 85.125, .25, Inds::lon));
    MatrixGeographicMap<double> mask_new(DividedRange::withEnds(29.625, 33.875, .25, Inds::lat),
                                         DividedRange::withEnds(74.875, 85.125, .25, Inds::lon));

    for (unsigned rr_new = 0; rr_new < mask_new.getLatitudes().count(); rr_new++)
      for (unsigned cc_new = 0; cc_new < mask_new.getLongitudes().count(); cc_new++) {
        unsigned numer = 0, denom = 0;
        double lat_total = 0, lon_total = 0; // for calculating exit direction
        // loop through elements within this frame
        for (Measure lat = mask_new.getLatitudes().getCellMin(rr_new) + mask.getLatitudes().getWidths() / 2;
             lat < mask_new.getLatitudes().getCellMax(rr_new); lat += mask.getLatitudes().getWidths())
          for (Measure lon = mask_new.getLongitudes().getCellMin(cc_new) + mask.getLongitudes().getWidths() / 2;
               lon < mask_new.getLongitudes().getCellMax(cc_new); lon += mask.getLongitudes().getWidths()) {
            denom++;

            if (lat < mask.getLatitudes().getMin() || lat >= mask.getLatitudes().getMax() ||
                lon < mask.getLongitudes().getMin() || lon >= mask.getLongitudes().getMax())
              continue;

            if (mask.getDouble(lat, lon) > 0) {
              numer++;

              // Follow this out of the cell
              queue< pair< pair<Measure, Measure>, double> > pending;
              pending.push(pair< pair<Measure, Measure>, double>(pair<Measure, Measure>(lat, lon), 1.0));
              
              while (!pending.empty()) {
                pair< pair<Measure, Measure>, double> llp = pending.front();
                pair<Measure, Measure> ll = llp.first;
                if (ll.first < mask_new.getLatitudes().getCellMin(rr_new) ||
                    ll.first >= mask_new.getLatitudes().getCellMax(rr_new) ||
                    ll.second < mask_new.getLongitudes().getCellMin(cc_new) ||
                    ll.second >= mask_new.getLongitudes().getCellMax(cc_new)) {
                  double lat_one = (ll.first - mask_new.getLatitudes().getCellCenter(rr_new)).getValue();
                  double lon_one = (ll.second - mask_new.getLongitudes().getCellCenter(cc_new)).getValue();
                  double length = sqrt(lat_one*lat_one + lon_one*lon_one);
                  lat_total += llp.second * lat_one / length;
                  lon_total += llp.second * lon_one / length;
                } else {
                  Measure lat0(Inds::lat), lon0(Inds::lon), lat1(Inds::lat), lon1(Inds::lon);
                  double portion0;
                  direction.getDirections(ll.first, ll.second, lat0, lon0, lat1, lon1, portion0);
                  if (portion0 > 0)
                    pending.push(pair< pair<Measure, Measure>, double>(pair<Measure, Measure>(lat0, lon0), portion0 * llp.second));
                  if (portion0 < 1)
                    pending.push(pair< pair<Measure, Measure>, double>(pair<Measure, Measure>(lat1, lon1), (1 - portion0) * llp.second));
                }
                    
                pending.pop();
              }
            }
          }
        
        mask_new.getCell(rr_new, cc_new) = numer / (double) denom;
        double angle = atan2(-lat_total, lon_total);
        if (angle < 0)
          angle += 2 * M_PI;
        direction_new.getCell(rr_new, cc_new) = angle;
      }

    DInfinityMap direction_new_map(direction_new);
    for (unsigned rr_new = 1; rr_new < mask_new.getLatitudes().count() - 1; rr_new++)
      for (unsigned cc_new = 1; cc_new < mask_new.getLongitudes().count() - 1; cc_new++) {
        if (mask_new.getCellConst(rr_new, cc_new) > 0) {
          // Do I flow into any that flow into me?
          // TODO: if I find these, have on hand the direction map made from the downsampled topography, and replace both with that, then run again
          checkBidirectional(rr_new, cc_new, rr_new - 1, cc_new - 1, direction_new_map);
          checkBidirectional(rr_new, cc_new, rr_new - 1, cc_new, direction_new_map);
          checkBidirectional(rr_new, cc_new, rr_new - 1, cc_new + 1, direction_new_map);
          checkBidirectional(rr_new, cc_new, rr_new, cc_new - 1, direction_new_map);
          checkBidirectional(rr_new, cc_new, rr_new, cc_new + 1, direction_new_map);
          checkBidirectional(rr_new, cc_new, rr_new + 1, cc_new - 1, direction_new_map);
          checkBidirectional(rr_new, cc_new, rr_new + 1, cc_new, direction_new_map);
          checkBidirectional(rr_new, cc_new, rr_new + 1, cc_new + 1, direction_new_map);
        }
      }

    direction_new.getValues().saveDelimited(argv[4], FileFormatter<double>::formatSimple, '\t');
    mask_new.getValues().saveDelimited(argv[3], FileFormatter<double>::formatSimple, '\t');
  } MPI_Finalize();
}

void checkBidirectional(unsigned rr0, unsigned cc0, unsigned rr1, unsigned cc1, DInfinityMap& direction_new) {
  double amount = direction_new.flowsIntoAmount(rr0, cc0, rr1, cc1) * direction_new.flowsIntoAmount(rr1, cc1, rr0, cc0);
  if (amount > .1)
    cout << "Bidirectional at " << rr0 << ", " << cc0 << " and " << rr1 << ", " << cc1 << ": " << amount << endl;
}
