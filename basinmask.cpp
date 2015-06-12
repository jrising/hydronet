#include <tools/hydro/DInfinityMap.h>

using namespace openworld;

void checkCell(unsigned rr_from, unsigned cc_from, unsigned rr_to, unsigned cc_to, DInfinityMap& map, queue< pair<unsigned, unsigned> >& pending, GeographicMap<bool>& mask);

// call as basinmask ang.tiff mask.tsv OUTLAT OUTLON LAT0 LAT1 DLAT LON0 LON1 DLON
int main(int argc, const char* argv[]) {
  MPI_Init(NULL,NULL); {
    int rank,size;
    MPI_Comm_rank(MCW,&rank);
    MPI_Comm_size(MCW,&size);

    double lat0 = atof(argv[5]);
    double lat1 = atof(argv[6]);
    double dlat = atof(argv[7]);
    double lon0 = atof(argv[8]);
    double lon1 = atof(argv[9]);
    double dlon = atof(argv[10]);
    
    DInfinityMap map(*MatrixGeographicMap<float>::loadTIFF(DividedRange::withEnds(lat0, lat1, dlat, Inds::lat),
                                                           DividedRange::withEnds(lon0, lon1, dlon, Inds::lon),
                                                           argv[1]), false);
    MatrixGeographicMap<bool> mask(map.getLatitudes(), map.getLongitudes());

    unsigned rrinit = map.getLatitudes().inRange(atof(argv[3])), ccinit = map.getLongitudes().inRange(atof(argv[4]));
    cout << "Flow to " << rrinit << ", " << ccinit << endl;

    mask.getCell(rrinit, ccinit) = true;

    queue< pair<unsigned, unsigned> > pending;
    pending.push(pair<unsigned, unsigned>(rrinit, ccinit));

    while (!pending.empty()) {
      pair<unsigned, unsigned> rc = pending.front();

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
