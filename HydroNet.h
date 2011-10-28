#ifndef HYDRONET_H
#define HYDRONET_H

namespace openworld {
  class HydroNode {
  protected:
    // structural
    list< pair<double, HydroNode*> > edges;

    double cellArea;
    double distAlong;
    double distAcross;
    double slope;
    double minToFlow;
    double maxVelocity;

    // time-varying state
    double volume;
    double conf;
    double volumeAfter;
    double confAfter;
    double confWeightAfter;
    bool stepped;

    double recursiveCalcResult;
    bool recursiveCalcDone;

    HydroNode(double cellArea, double distAlong, double distAcross, double slope, double minToFlow, double maxVelocity) {
      this->cellArea = cellArea;
      this->distAlong = distAlong;
      this->distAcross = distAcross;
      this->slope = slope;
      this->minToFlow = minToFlow;
      this->maxVelocity = maxVelocity;

      volume = conf = volumeAfter = confAfter = confWeightAfter = 0;
      stepped = false;
    }

  public:
    // ignore keys of map
    void setEdges(map<unsigned, pair<double, HydroNode*>> edges, double divide = 1.0) {
      this->edges.clear();

      double total = 0;
      map<unsigned, pair<double, HydroNode*>>::iterator it;
      for (it = edges.begin(); it != edges.end(); it++) {
        edges.push_back(pair<double, HydroNode*>(it->second.first / divide, it->second.second));
        total += it->second.first / divide;
      }

      if (total != 1.0)
        throw runtime_error("Edges do not sum to 1");
    }

    virtual double calcVelocity() = 0;

    double addStepVolume(double volume, double conf, double confWeight) {
      volumeAfter += volume;
      confAfter += conf * confWeight;
      confWeightAfter += confWeight;
    }

    virtual void stepModel(double dt) {
      if (stepped)
        return;

      if (volume >= minToFlow) {
        double portion = calcVelocity() * step / distAlong;

        list< pair<double, HydroNode*> >::iterator it;
        for (it = edges.begin(); it != edges.end(); it++)
          it->second->addStepVolume(portion * it->first * volume, portion * it->first * conf, portion * it->first);

        addStepVolume((1 - portion) * volume, (1 - portion) * it->first * conf, (1 - portion) * it->first);
      } else
        addStepVolume(volume, conf, 1);

      stepped = true;

      // Recurse on all edges
      list< pair<double, HydroNode*> >::iterator it;
      for (it = edges.begin(); it != edges.end(); it++)
        it->second->stepModel(dt);
    }

    virtual void stepPost() {
      if (!stepped)
        return;

      volume = volumeAfter;
      if (confWeightAfter == 0)
        conf = 0;
      else
        conf = confAfter / confWeightAfter;

      volumeAfter = confAfter = confWeightAfter = 0;
      stepped = false;

      // Recurse on all edges
      list< pair<double, HydroNode*> >::iterator it;
      for (it = edges.begin(); it != edges.end(); it++)
        it->second->stepPost();
    }

    // Slope is drop / run
    static double calcManning(double coeff, double radius, double slope) {
      if (radius <= 0 || slope <= 0)
        return 0.0;
      double vel = (1 / coeff) * pow(radius, 2.0/3.0) * sqrt(slope);

      return vel;
    }

    // recursive calculations
    void resetRecursiveCalculation() {
      if (!recursiveCalcDone)
        return;

      recursiveCalcDone = false;

      list< pair<double, HydroNode*> >::iterator it;
      for (it = edges.begin(); it != edges.end(); it++)
        if (it->second->recursiveCalcDone)
          it->second->resetRecursiveCalculation();
    }

    double recursiveMaximumStep() {
      recursiveCalcDone = true;

      double vel = calcVelocity();
      if (vel > 0)
        recurisveCalcResult = distAlong / vel;
      else
        recursiveCalcResult = numerical_limits<double>::max();

      list< pair<double, HydroNode*> >::iterator it;
      for (it = edges.begin(); it != edges.end(); it++)
        if (!it->second->recursiveCalcDone) { // otherwise, ignore
          double step = it->second->recursiveMaximumStep();
          if (step < recursiveCalcResult)
            recursiveCalcResult = step;
        }

      return recurisveCalcResult;
    }
  }; 

  class HydroSurfaceNode {
  private:
    double minSurfaceToFlow = .1;
    double maxSurfaceVelocity = 1;

  public:
    HydroSurfaceNode(double cellArea, double distAlong, double distAcross, double slope)
      : HydroNode(cellArea, distAlong, distAcross, slope, minSurfaceToFlow, maxSurfaceVelocity) {
    }

    virtual double calcVelocity() {
      return calcManningSurface(volume);
    }

    static double calcManningSurface(double volume) {
      if (volume < minToFlow)
        return 0;
      
      double height = volume / cellArea;
      double vel = calcManning(.025, height, max(height / (2 * distAlong), slope));
      return min(vel, maxVelocity); //90 * sqrt(height / 2)); //(water terminal velocity)
    }
  };

  class HydroRiverNode {
  private:
    double minRiverToFlow = .001;
    double maxRiverVelocity = 5;

  public:
    HydroRiverNode(double cellArea, double distAlong, double distAcross, double slope)
      : HydroNode(cellArea, distAlong, distAcross, slope, minRiverToFlow, maxRiverVelocity) {
    }

    virtual double calcVelocity() {
      return calcManningRiver(volume);
    }

    static double calcRiverWidth(double volume) {
      return 2 * sqrt((3 / M_PI) * volume / distAlong);
    }

    static double calcManningRiver(double volume) {
      if (volume < minToFlow)
        return 0;

      // 1/3 of circle calculation
      double width = calcRiverWidth(volume);
      if (width <= 0)
        return 0;

      double r, height;
      if (width > distAcross) {
        r = volume / cellArea;
        height = r;
      } else {
        r = (volume / distAlong) / (M_PI * width / 3);
        height = width / 2;
      }

      double vel = calcManning(.033, r, max(height / (2 * distAlong), slope));
      return min(vel, maxVelocity); //90 * sqrt(height)); // terminal velocity
    }
  };

  class HydroOutputNode {
  public:
    HydroOutputNode()
      : HydroNode(0, 0, 0, 0, 0, 0) {
    }
    
    virtual void stepModel(double dt) {
      // ignore
    }

    virtual void stepPost() {
      // ignore
    }
  };
    
  class HydroNet {
  protected:
    GeographicMap<double> mask_coarse;
    map<unsigned, HydroNode*> surfaces_coarse;

  public:
    HydroNet(GeographicMap<double>& mask_coarse)
      : mask_coarse(mask_coarse) {
    }

    HydroOutputNode* generate(DInfinityMap& direction_fine, GeographicMap<bool>& mask_fine, GeographicMap<double>& slope_fine) {
      HydroOutputNode* out = new HydroOutputNode();

      map<unsigned, HydroNode*> rivers;

      surfaces_coarse.clear();

      // Loop over all surface cells
      for (unsigned rr = 0; rr < mask_coarse.getLatitudes().count(); rr++)
        for (unsigned cc = 0; cc < mask_coarse.getLongitudes().count(); cc++) {
          double total_dist = 0, total_slope = 0; weight_slope = 0;
          unsigned paths = 0;

          map<unsigned, pair<double, HydroNode*>> edges;
                    
          // Loop over all elevation cells within
          for (Measure lat = mask_coarse.getLatitudes().getCellMin(rr) + mask_fine.getLatitudes().getWidths() / 2;
               lat < mask_coarse.getLatitudes().getCellMax(rr); lat += mask_fine.getLatitudes().getWidths())
            for (Measure lon = mask_coarse.getLongitudes().getCellMin(cc) + mask_fine.getLongitudes().getWidths() / 2;
                 lon < mask_coarse.getLongitudes().getCellMax(cc); lon += mask_fine.getLongitudes().getWidths()) {
              
              if (lat < mask_fine.getLatitudes().getMin() || lat >= mask_fine.getLatitudes().getMax() ||
                  lon < mask_fine.getLongitudes().getMin() || lon >= mask_fine.getLongitudes().getMax() ||
                  mask_fine.getDouble(lat, lon) == 0)
                continue;

              paths++;

              // Follow this out of the cell
              queue< pair< pair<Measure, Measure>, double> > pending;
              pending.push(pair< pair<Measure, Measure>, double>(pair<Measure, Measure>(lat, lon), 1.0));
              
              while (!pending.empty()) {
                pair< pair<Measure, Measure>, double> llp = pending.front();
                pair<Measure, Measure> ll = llp.first;
                if (ll.first < mask_coarse.getLatitudes().getCellMin(rr) ||
                    ll.first >= mask_coarse.getLatitudes().getCellMax(rr) ||
                    ll.second < mask_coarse.getLongitudes().getCellMin(cc) ||
                    ll.second >= mask_coarse.getLongitudes().getCellMax(cc)) {
                  // is this the output node?
                  if (mask_fine.getDouble(ll.first, ll.second) == 0)
                    rivers[mask_fine.getIndex(ll.first, ll.second)] = out;

                  // connect this to the river out
                  if (rivers.find(mask_fine.getIndex(ll.first, ll.second)) == map::end)
                    generateRiver(ll.first, ll.second, rivers, direction_fine, mask_fine, slope_fine, out);
                  
                  if (edges.find(mask_fine.getIndex(ll.first, ll.second)) == map::end)
                    edges[mask_fine.getIndex(ll.first, ll.second)] = 
                      pair<double, HydroNode*>(llp.second, rivers[mask_fine.getIndex(ll.first, ll.second)]);
                  else
                    edges[mask_fine.getIndex(ll.first, ll.second)].first += llp.second;
                } else {
                  total_dist += llp.second * calcDistance(rr - 1, cc - 1, rr + 1, cc + 1) / (3 * M_SQRT2);
                  total_slope += llp.second * slope_fine.getDouble(lat, lon);
                  weight_slope += llp.second;

                  Measure lat0(Inds::lat), lon0(Inds::lon), lat1(Inds::lat), lon1(Inds::lon);
                  double portion0;

                  direction_fine.getDirections(ll.first, ll.second, lat0, lon0, lat1, lon1, portion0);
                  if (portion0 > 0)
                    pending.push(pair< pair<Measure, Measure>, double>(pair<Measure, Measure>(lat0, lon0), portion0 * llp.second));
                  if (portion0 < 1)
                    pending.push(pair< pair<Measure, Measure>, double>(pair<Measure, Measure>(lat1, lon1), (1 - portion0) * llp.second));
                }

                pending.pop();
              }
            }

          HydroNode* surface = new HydroSurfaceNode(mask_coarse.calcArea(rr, cc) * mask_coarse.getCellConst(rr, cc),
                                                    total_dist / paths, calcDistance(rr - 1, cc - 1, rr + 1, cc + 1) / (3 * M_SQRT2),
                                                    total_slope / weight_slope, false, minSurfaceToFlow, maxSurfaceVelocity);

          // divide edges weights by #paths
          surface->setEdges(edges, paths);
          surfaces_coarse[mask_coarse.getIndex(rr, cc)] = surface;
        }

      return out;
    }
    
    void generateRiver(Measure lat, Measure lon, map<unsigned, HydroNode*> rivers, DInfinityMap& direction_fine, GeographicMap<bool>& mask_fine, GeographicMap<double>& slope_fine, HydroOutputNode* out) {
      map<unsigned, pair<double, HydroNode*>> edges;
      double total_dist = 0, total_slope = 0; weight_slope = 0;

      // Follow this out of the cell
      queue< pair< pair<Measure, Measure>, double> > pending;
      pending.push(pair< pair<Measure, Measure>, double>(pair<Measure, Measure>(lat, lon), 1.0));
              
      while (!pending.empty()) {
        pair< pair<Measure, Measure>, double> llp = pending.front();
        pair<Measure, Measure> ll = llp.first;
        if (ll.first < mask_coarse.getLatitudes().getCellMin(rr) ||
            ll.first >= mask_coarse.getLatitudes().getCellMax(rr) ||
            ll.second < mask_coarse.getLongitudes().getCellMin(cc) ||
            ll.second >= mask_coarse.getLongitudes().getCellMax(cc)) {
          // is this the output node?
          if (mask_fine.getDouble(ll.first, ll.second) == 0)
            rivers[mask_fine.getIndex(ll.first, ll.second)] = out;

          if (rivers.find(mask_fine.getIndex(ll.first, ll.second)) == map::end)
            generateRiver(ll.first, ll.second, rivers, direction_fine, mask_fine, slope_fine, out);
          
          if (edges.find(mask_fine.getIndex(ll.first, ll.second)) == map::end)
            edges[mask_fine.getIndex(ll.first, ll.second)] = 
              pair<double, HydroNode*>(llp.second, rivers[mask_fine.getIndex(ll.first, ll.second)]);
          else
            edges[mask_fine.getIndex(ll.first, ll.second)].first += llp.second;
        } else {
          total_dist += llp.second * calcDistance(rr - 1, cc - 1, rr + 1, cc + 1) / (3 * M_SQRT2);
          total_slope += llp.second * slope_fine.getDouble(lat, lon);
          weight_slope += llp.second;

          Measure lat0(Inds::lat), lon0(Inds::lon), lat1(Inds::lat), lon1(Inds::lon);
          double portion0;

          direction_fine.getDirections(ll.first, ll.second, lat0, lon0, lat1, lon1, portion0);
          if (portion0 > 0)
            pending.push(pair< pair<Measure, Measure>, double>(pair<Measure, Measure>(lat0, lon0), portion0 * llp.second));
          if (portion0 < 1)
            pending.push(pair< pair<Measure, Measure>, double>(pair<Measure, Measure>(lat1, lon1), (1 - portion0) * llp.second));
        }

        pending.pop();
      }

      // actually make the river and add it to rivers
      HydroNode* river = new HydroRiverNode(mask_coarse.calcArea(rr, cc) * mask_coarse.getCellConst(rr, cc),
                                            total_dist / paths, calcDistance(rr - 1, cc - 1, rr + 1, cc + 1) / (3 * M_SQRT2),
                                            total_slope / weight_slope, true, minRiverToFlow, maxRiverVelocity);
      
      river->setEdges(edges);
      rivers[mask_fine.getIndex(rr, cc)] = river;
    }

    double calculateMaximumStep() {
      double maxstep = numerical_limits<double>::max();

      map<unsigned, HydroNode*>::iterator it;
      for (it = surfaces_coarse.begin(); it != surfaces_coarse.end(); it++)
        it->second->resetRecursiveCalculation();

      for (it = surfaces_coarse.begin(); it != surfaces_coarse.end(); it++) {
        double step = it->second->recursiveMaximumStep();
        if (step < maxstep)
          maxstep = step;
      }
      
      return maxstep;
    }

    void step(double dt, GeographicMap<double>& changes, GeographicMap<double>& changesConf, double scale) {
      map<unsigned, HydroNode*>::iterator it;
      for (it = surfaces_coarse.begin(); it != surfaces_coarse.end(); it++)
        it->second->stepModel(dt);

      for (unsigned rr = 0; rr < mask_coarse.getLatitudes().count(); rr++)
        for (unsigned cc = 0; cc < mask_coarse.getLongitudes().count(); cc++) {
          double fraction = mask_coarse.getCellConst(rr, cc);
          if (fraction == 0)
            continue;

          surfaces_coarse[mask_coarse.getIndex(rr, cc)].addStepVolume(changes.getCell(rr, cc) * fraction * scale, changesConf.getCell(rr, cc), 1.0);
        }

      for (it = surfaces_coarse.begin(); it != surfaces_coarse.end(); it++)
        it->second->stepPost();
    }
  };
}

#endif
