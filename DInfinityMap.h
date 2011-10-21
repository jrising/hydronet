#ifndef DINFINITY_MAP_H
#define DINFINITY_MAP_H

#include <datastr/GeographicMap.h>

namespace openworld {
  class DInfinityMap : public GeographicMap<double> {
  protected:
    int latsucc;

  public:
    template<class S>
    DInfinityMap(const GeographicMap<S>& src, bool latincdown = false)
      : GeographicMap<double>(src) {

      if (latincdown)
        latsucc = -1;
      else
        latsucc = 1;
    }

    void getDirections(unsigned rr, unsigned cc, unsigned& rr0, unsigned& cc0, unsigned& rr1, unsigned& cc1, double& portion0) const {
      // rr0, cc0 is always straight, rr1, cc1 is diagonal
      // portion1 = 1 - portion0

      double dir = getCellConst(rr, cc);
      if (dir < 0 || dir > 2 * PI)
        dir = 2 * PI * ((double) rand() / (double) RAND_MAX);
          
      if (dir < M_PI / 4) {
        rr0 = rr;
        cc0 = cc1 = cc + 1;
        rr1 = rr - latsucc;
        portion0 = 1 - dir / (M_PI / 4);
      } else if (dir < M_PI / 2) {
        rr0 = rr1 = rr - latsucc;
        cc0 = cc;
        cc1 = cc + 1;
        portion0 = (dir - M_PI / 4) / (M_PI / 4);
      } else if (dir < 3 * M_PI / 4) {
        rr0 = rr1 = rr - latsucc;
        cc0 = cc;
        cc1 = cc - 1;
        portion0 = 1 - (dir - M_PI / 2) / (M_PI / 4);
      } else if (dir < M_PI) {
        rr0 = rr;
        cc0 = cc1 = cc - 1;
        rr1 = rr - latsucc;
        portion0 = (dir - 3 * M_PI / 4) / (M_PI / 4);
      } else if (dir < 5 * M_PI / 4) {
        rr0 = rr;
        cc0 = cc1 = cc - 1;
        rr1 = rr + latsucc;
        portion0 = 1 - (dir - M_PI) / (M_PI / 4);
      } else if (dir < 3 * M_PI / 2) {
        rr0 = rr1 = rr + latsucc;
        cc0 = cc;
        cc1 = cc - 1;
        portion0 = (dir - 5 * M_PI / 4) / (M_PI / 4);
      } else if (dir < 7 * M_PI / 4) {
        rr0 = rr1 = rr + latsucc;
        cc0 = cc;
        cc1 = cc + 1;
        portion0 = 1 - (dir - 3 * M_PI / 2) / (M_PI / 4);
      } else {
        rr0 = rr;
        cc0 = cc1 = cc + 1;
        rr1 = rr + latsucc;
        portion0 = (dir - 7 * M_PI / 4) / (M_PI / 4);
      }

      if (portion0 < 1e-6)
        portion0 = 0;
      if (1 - portion0 < 1e-6)
        portion0 = 1;
    }

    bool flowsInto(unsigned rr_from, unsigned cc_from, unsigned rr_to, unsigned cc_to) {
      unsigned rr0, cc0, rr1, cc1;
      double portion0;

      getDirections(rr_from, cc_from, rr0, cc0, rr1, cc1, portion0);
      return ((rr0 == rr_to && cc0 == cc_to && portion0 > 0) || (rr1 == rr_to && cc1 == cc_to && portion0 < 1));
    }

    double flowsIntoAmount(unsigned rr_from, unsigned cc_from, unsigned rr_to, unsigned cc_to) {
      unsigned rr0, cc0, rr1, cc1;
      double portion0;

      getDirections(rr_from, cc_from, rr0, cc0, rr1, cc1, portion0);
      if (rr0 == rr_to && cc0 == cc_to)
        return portion0;
      if (rr1 == rr_to && cc1 == cc_to)
        return 1 - portion0;
      return 0;
    }

    void getDirections(Measure lat, Measure lon, Measure& lat0, Measure& lon0, Measure& lat1, Measure& lon1, double& portion0) const {
      unsigned rr0, cc0, rr1, cc1;
      getDirections(getLatitudes().inRange(lat), getLongitudes().inRange(lon), rr0, cc0, rr1, cc1, portion0);
      
      lat0 = getLatitudes().getCellCenter(rr0);
      lon0 = getLongitudes().getCellCenter(cc0);
      lat1 = getLatitudes().getCellCenter(rr1);
      lon1 = getLongitudes().getCellCenter(cc1);
    }
  };
}

#endif
