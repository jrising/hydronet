#include "DInfinityMap.h"

using namespace openworld;

void checkCell(unsigned rr_from, unsigned cc_from, unsigned rr_to, unsigned cc_to, DInfinityMap& map, queue< pair<unsigned, unsigned> >& pending, GeographicMap<bool>& mask);

// call as basinmask ang.tiff mask.tsv 32.2341 -2.324
int main(int argc, const char* argv[]) {
  MPI_Init(NULL,NULL); {
	int rank,size;
	MPI_Comm_rank(MCW,&rank);
	MPI_Comm_size(MCW,&size);

    DInfinityMap map(*MatrixGeographicMap<float>::loadTIFF(DividedRange::withEnds(29.74583, 33.9125, 0.008333334, Inds::lat),
                                                           DividedRange::withEnds(74.77084, 85.17084, 0.008333334, Inds::lon),
                                                           argv[1]), false);
    MatrixGeographicMap<bool> mask(map.getLatitudes(), map.getLongitudes());

    unsigned rrinit = map.getLatitudes().inRange(atof(argv[3])), ccinit = map.getLongitudes().inRange(atof(argv[4]));
    cout << "Flow to " << rrinit << ", " << ccinit << endl;

    mask.getCell(rrinit, ccinit) = true;

    queue< pair<unsigned, unsigned> > pending;
    pending.push(pair<unsigned, unsigned>(rrinit, ccinit));

    //unsigned furthest_cc = 0, furthest_rr = 0;

    while (!pending.empty()) {
      pair<unsigned, unsigned> rc = pending.front();

      //cout << rc.first << ", " << rc.second << ": " << map.getCellConst(rc.first, rc.second) << endl;

      /*if (rc.second > furthest_cc) {
        furthest_rr = rc.first;
        furthest_cc = rc.second;
        }*/

      checkCell(rc.first - 1, rc.second - 1, rc.first, rc.second, map, pending, mask);
      checkCell(rc.first - 1, rc.second, rc.first, rc.second, map, pending, mask);
      checkCell(rc.first - 1, rc.second + 1, rc.first, rc.second, map, pending, mask);
      checkCell(rc.first, rc.second - 1, rc.first, rc.second, map, pending, mask);
      checkCell(rc.first, rc.second + 1, rc.first, rc.second, map, pending, mask);
      checkCell(rc.first + 1, rc.second - 1, rc.first, rc.second, map, pending, mask);
      checkCell(rc.first + 1, rc.second, rc.first, rc.second, map, pending, mask);
      checkCell(rc.first + 1, rc.second + 1, rc.first, rc.second, map, pending, mask);

      pending.pop();
    }

    /*for (unsigned rr = furthest_rr - 2; rr <= furthest_rr + 2; rr++) {
      for (unsigned cc = furthest_cc - 2; cc <= furthest_cc + 2; cc++)
        cout << map.getCellConst(rr, cc) << ":" << mask.getCellConst(rr, cc) << " ";
      cout << endl;
      }

    // figure out where it would go!
    unsigned rr = furthest_rr, cc = furthest_cc + 1;
    while (rr >= 0 || rr < map.getLatitudes().count() || cc > 0 || cc < map.getLongitudes().count()) {
      cout << rr << ", " << cc << ": " << map.getCellConst(rr, cc) << endl;
      mask.getCell(rr, cc) = -1;

      unsigned rr0, cc0, rr1, cc1;
      double portion0;
      map.getDirections(rr, cc, rr0, cc0, rr1, cc1, portion0);
      cout << "  " << rr0 << ", " << cc0 << " and " << rr1 << ", " << cc1 << " - " << portion0 << endl;
      if ((mask.getCellConst(rr0, cc0) == 1 && portion0 > 0) || (mask.getCellConst(rr1, cc1) == 1 && portion0 < 1)) {
        cout << "Flows in!" << endl;
        break;
      } else if ((portion0 == 0 || mask.getCellConst(rr0, cc0) == -1) && (portion0 == 1 || mask.getCellConst(rr1, cc1) == -1)) {
        cout << "Flows in circles!" << endl;
        break;
      } else if ((mask.getCellConst(rr0, cc0) == -1 && portion0 < 1) || portion0 == 0) {
        rr = rr1;
        cc = cc1;
      } else {
        rr = rr0;
        cc = cc0;
      }
      }*/

    mask.getValues().saveDelimited(argv[2], FileFormatter<bool>::formatSimple, '\t');
  } MPI_Finalize();
}

void checkCell(unsigned rr_from, unsigned cc_from, unsigned rr_to, unsigned cc_to, DInfinityMap& map, queue< pair<unsigned, unsigned> >& pending, GeographicMap<bool>& mask) {
  if (rr_from < 0 || rr_from >= map.getLatitudes().count() || cc_from < 0 || cc_from >= map.getLongitudes().count())
    return;
  if (mask.getCellConst(rr_from, cc_from))
    return;
  if (map.flowsInto(rr_from, cc_from, rr_to, cc_to)) {
    mask.getCell(rr_from, cc_from) = true;
    pending.push(pair<unsigned, unsigned>(rr_from, cc_from));
  }
}
