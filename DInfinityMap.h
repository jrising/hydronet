#ifndef DINFINITY_MAP_H
#define DINFINITY_MAP_H

#include <datastr/GeographicMap.h>

namespace openworld {
  class DInfinityMap : public GeographicMap<double> {
  public:
    template<class S>
    DInfinityMap(const GeographicMap<S>& src)
      : GeographicMap<double>(src) {
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
        rr1 = rr - 1;
        portion0 = 1 - dir / (M_PI / 4);
      } else if (dir < M_PI / 2) {
        rr0 = rr1 = rr - 1;
        cc0 = cc;
        cc1 = cc + 1;
        portion0 = (dir - M_PI / 4) / (M_PI / 4);
      } else if (dir < 3 * M_PI / 4) {
        rr0 = rr1 = rr - 1;
        cc0 = cc;
        cc1 = cc - 1;
        portion0 = 1 - (dir - M_PI / 2) / (M_PI / 4);
      } else if (dir < M_PI) {
        rr0 = rr;
        cc0 = cc1 = cc - 1;
        rr1 = rr - 1;
        portion0 = (dir - 3 * M_PI / 4) / (M_PI / 4);
      } else if (dir < 5 * M_PI / 4) {
        rr0 = rr;
        cc0 = cc1 = cc - 1;
        rr1 = rr + 1;
        portion0 = 1 - (dir - M_PI) / (M_PI / 4);
      } else if (dir < 3 * M_PI / 2) {
        rr0 = rr1 = rr + 1;
        cc0 = cc;
        cc1 = cc - 1;
        portion0 = (dir - 5 * M_PI / 4) / (M_PI / 4);
      } else if (dir < 7 * M_PI / 4) {
        rr0 = rr1 = rr + 1;
        cc0 = cc;
        cc1 = cc + 1;
        portion0 = 1 - (dir - 3 * M_PI / 2) / (M_PI / 4);
      } else {
        rr0 = rr;
        cc0 = cc1 = cc + 1;
        rr1 = rr + 1;
        portion0 = (dir - 7 * M_PI / 4) / (M_PI / 4);
      }
    }

    bool flowsInto(unsigned rr_from, unsigned cc_from, unsigned rr_to, unsigned cc_to) {
      unsigned rr0, cc0, rr1, cc1;
      double portion0;

      map.getDirections(rr_from, cc_from, rr0, cc0, rr1, cc1, portion0);
      return ((rr0 == rr_to && cc0 == cc_to) || (rr1 == rr_to && cc1 == cc_to))
    }
  };
}

#endif
