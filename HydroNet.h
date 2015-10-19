#ifndef HYDRONET_H
#define HYDRONET_H

#include <float.h>
#include <utils/ToString.h>
#include <memory/SerializationTools.h>
#include <list>
#include <queue>
#include <measure/Measure.h>
#include <datastr/GeographicMap.h>
#include <tools/hydro/DInfinityMap.h>

using namespace std;

namespace openworld {
  void* nodeListConstructor(istream& in, PointerReference& reference);

  class HydroNode : public IPointerSerializable {
  protected:
    // structural
    list< pair<double, HydroNode*> > edges;

    double nodeArea;
    double distAlong;
    double distAcross;
    double slope;
    double minToFlow;
    double maxVelocity;

    // meta-data
    Measure latitude;
    Measure longitude;

    // time-varying state
    double precipVolume;
    double meltVolume;
    double conf;
    double precipVolumeAfter;
    double meltVolumeAfter;
    double confAfter;
    double confWeightAfter;
    bool stepped;

    double evaporationRate;

    double recursiveCalcResult;
    bool recursiveCalcDone;

    HydroNode(double nodeArea, double distAlong, double distAcross, double slope, double minToFlow, double maxVelocity, Measure latitude, Measure longitude)
      : latitude(latitude), longitude(longitude) {
      this->nodeArea = nodeArea;
      this->distAlong = distAlong;
      this->distAcross = distAcross;
      this->slope = slope;
      this->minToFlow = minToFlow;
      this->maxVelocity = maxVelocity;
      this->evaporationRate = 0;

      precipVolume = meltVolume = conf = precipVolumeAfter = meltVolumeAfter = confAfter = confWeightAfter = 0;
      stepped = false;
    }

  public:
    virtual HydroNode* clone(map<HydroNode*, HydroNode*>& translate) = 0;

    Measure getLatitude() {
      return latitude;
    }

    Measure getLongitude() {
      return longitude;
    }

    double getPrecipVolume() {
      return precipVolume;
    }

    double getMeltVolume() {
      return meltVolume;
    }

    // ignore keys of map
    void setEdges(map<unsigned, pair<double, HydroNode*> >& edges, double divide = 1.0) {
      this->edges.clear();
      if (divide == 0)
        return;

      double total = 0;
      map<unsigned, pair<double, HydroNode*> >::iterator it;
      for (it = edges.begin(); it != edges.end(); it++) {
        this->edges.push_back(pair<double, HydroNode*>(it->second.first / divide, it->second.second));
        total += it->second.first / divide;
      }

      if (total < 1.0 - FLT_EPSILON * edges.size() || total > 1.0 + FLT_EPSILON * edges.size())
        throw runtime_error("Edges do not sum to 1: " + ToString::flong(total) + " for " + ToString::base10(divide));
    }

    virtual double calcVelocity() = 0;
    virtual double calcEvaporationRate(GeographicMap<double>& temps, GeographicMap<double>& snowCover, GeographicMap<double>& relhums, GeographicMap<double>& lights, double fracyear, double surfaceEvaporationFactor, double riverEvaporationFactor) = 0;

    void addStepVolume(double precipVolume, double meltVolume, double conf, double confWeight) {
      precipVolumeAfter += precipVolume;
      meltVolumeAfter += meltVolume;
      confAfter += conf * confWeight;
      confWeightAfter += confWeight;
    }

    virtual void stepModel(double dt) {
      if (stepped)
        return;

      double afterEvapPrecipVolume = precipVolume, afterEvapMeltVolume = meltVolume;
      if (evaporationRate > 0) {
        double evapVolume = (0.001 * evaporationRate * nodeArea * dt) / DividedRange::toTimespan(1).getValue();

        afterEvapPrecipVolume = max(0.0, precipVolume - precipVolume * evapVolume / (precipVolume + meltVolume));
        afterEvapMeltVolume = max(0.0, meltVolume - meltVolume * evapVolume / (precipVolume + meltVolume));
      }

      if (afterEvapPrecipVolume + afterEvapMeltVolume >= minToFlow) {
        double portion = calcVelocity() * dt / distAlong;

        list< pair<double, HydroNode*> >::iterator it;
        for (it = edges.begin(); it != edges.end(); it++)
          it->second->addStepVolume(portion * it->first * afterEvapPrecipVolume, portion * it->first * afterEvapMeltVolume, portion * it->first * conf, portion * it->first);

        addStepVolume((1 - portion) * afterEvapPrecipVolume, (1 - portion) * afterEvapMeltVolume, (1 - portion) * it->first * conf, (1 - portion) * it->first);
      } else
        addStepVolume(afterEvapPrecipVolume, afterEvapMeltVolume, conf, 1);

      stepped = true;

      // Recurse on all edges
      list< pair<double, HydroNode*> >::iterator it;
      for (it = edges.begin(); it != edges.end(); it++)
        it->second->stepModel(dt);
    }

    virtual void stepPost() {
      if (!stepped)
        return;

      precipVolume = precipVolumeAfter;
      meltVolume = meltVolumeAfter;
      if (confWeightAfter == 0)
        conf = 0;
      else
        conf = confAfter / confWeightAfter;

      precipVolumeAfter = meltVolumeAfter = confAfter = confWeightAfter = 0;
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

    // From "Simpliﬁed versions for the Penman evaporation equation using routine weather data" by Valiantzas, 2006
    // temp in C, relhum between 0 and 1, light between 0 and 1 (portion of daylight hours light), lat in degrees, month = 1 for jan
    // the result is mm/day
    static double calcEvaporation(double temp, double relhum, double light, Measure lat, int month) {
      return 11.63 * exp(0.077 * temp) / 30;
    }

    /*static double calcEvaporation(double temp, double relhum, double light, Measure lat, int month) {
      double latrad = lat.getValue() * M_PI / 180;
      double N = 4 * latrad * sin(.053 * month - 1.65) + 12; // 12 at equator
      double R_A;
      //if (abs(latrad) > 23.5 * M_PI / 180)
      R_A = 3 * N * sin(.131*N - .95*latrad); // 36 at equator
      //else
      //R_A = 118 * pow(N, .2) * sin(.131*N - .2*latrad); // 193.9626 at equator -- seems way too high!

      double R_S = R_A * (.5 + .25 * light); // 27 at full bright

      // TODO: add 0:00012 * elevation

      return max(0.0, .047*R_S*sqrt(temp + 9.5) - 2.4*pow(R_S/R_A, 2) + .09*(temp + 20)*(1 - relhum));
      // 9.1424342 at no humid, 20 deg C
      // 7.3424342 at 50% humid, 20 deg C
      // 3.4613207 at 50% humid, 0 deg C
      // 2.0075471 at 100% humid (and no full-bright hours), 0 deg C
      }*/

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
        recursiveCalcResult = distAlong / vel;
      else
        recursiveCalcResult = FLT_MAX;

      list< pair<double, HydroNode*> >::iterator it;
      for (it = edges.begin(); it != edges.end(); it++)
        if (!it->second->recursiveCalcDone) { // otherwise, ignore
          double step = it->second->recursiveMaximumStep();
          if (step < recursiveCalcResult)
            recursiveCalcResult = step;
        }

      return recursiveCalcResult;
    }

    list<pair<pair<HydroNode*, HydroNode*>, double> > recursiveGetAllEdges() {
      recursiveCalcDone = true;
      list<pair<pair<HydroNode*, HydroNode*>, double> > result;

      list< pair<double, HydroNode*> >::iterator it;
      for (it = edges.begin(); it != edges.end(); it++) {
        result.push_back(pair<pair<HydroNode*, HydroNode*>, double>(pair<HydroNode*, HydroNode*>(this, it->second), it->first));

        if (!it->second->recursiveCalcDone) { // otherwise, don't recurse
          list<pair<pair<HydroNode*, HydroNode*>, double> > resultSet = it->second->recursiveGetAllEdges();
          result.insert(result.end(), resultSet.begin(), resultSet.end());
        }
      }

      return result;
    }

    void recursiveUpdateEvaporationRate(GeographicMap<double>& temps, GeographicMap<double>& snowCover, GeographicMap<double>& relhums, GeographicMap<double>& lights, double fracyear, double surfaceEvaporationFactor, double riverEvaporationFactor) {
      recursiveCalcDone = true;

      evaporationRate = calcEvaporationRate(temps, snowCover, relhums, lights, fracyear, surfaceEvaporationFactor, riverEvaporationFactor);

      list< pair<double, HydroNode*> >::iterator it;
      for (it = edges.begin(); it != edges.end(); it++)
        if (!it->second->recursiveCalcDone) { // otherwise, ignore
          it->second->recursiveUpdateEvaporationRate(temps, snowCover, relhums, lights, fracyear, surfaceEvaporationFactor, riverEvaporationFactor);
        }
    }

    void translateEdges(HydroNode& copy, map<HydroNode*, HydroNode*>& translate) {
      edges.clear();

      for (list< pair<double, HydroNode*> >::iterator it = copy.edges.begin(); it != copy.edges.end(); it++)
        edges.push_back(pair<double, HydroNode*>(it->first, HydroNode::getCopy(it->second, translate)));
    }

    static HydroNode* getCopy(HydroNode* node, map<HydroNode*, HydroNode*>& translate) {
      map<HydroNode*, HydroNode*>::iterator it = translate.find(node);
      if (it != translate.end())
        return it->second;

      HydroNode* copy = node->clone(translate);
      translate[node] = copy;

      return copy;
    }

    // Serializable protocol

    virtual ostream& streamInsert(ostream& os, PointerTracker& tracker) const {
      throw logic_error("HydroNode::streamInsert not implemented yet!");
    }

    virtual istream& streamExtract(istream& in, PointerReference& reference) {
      throw logic_error("HydroNode::streamExtract not implemented yet!");
    }

    static HydroNode* streamExtractPointer(istream& in, PointerReference& reference) {
      throw logic_error("HydroNode::streamExtractPointer not implemented yet!");
    }
  };

  class HydroSurfaceNode : public HydroNode {
  private:
    static const double minSurfaceToFlow = .1;
    static const double maxSurfaceVelocity = 1;

  public:
    HydroSurfaceNode(double nodeArea, double distAlong, double distAcross, double slope, Measure latitude, Measure longitude)
      : HydroNode(nodeArea, distAlong, distAcross, slope, minSurfaceToFlow, maxSurfaceVelocity, latitude, longitude) {
    }

    virtual HydroNode* clone(map<HydroNode*, HydroNode*>& translate) {
      HydroNode* copy = new HydroSurfaceNode(*this);
      copy->translateEdges(*this, translate);

      return copy;
    }

    virtual double calcVelocity() {
      return calcManningSurface(precipVolume + meltVolume);
    }

    virtual double calcEvaporationRate(GeographicMap<double>& temps, GeographicMap<double>& snowCover, GeographicMap<double>& relhums, GeographicMap<double>& lights, double fracyear, double surfaceEvaporationFactor, double riverEvaporationFactor) {
      double evapRate = HydroNode::calcEvaporation(temps.getDouble(latitude, longitude) - 273.15, relhums.getDouble(latitude, longitude), lights.getDouble(latitude, longitude), latitude, max(fracyear * 12 + 1, 12.0));
      return surfaceEvaporationFactor * evapRate * (1.0 - snowCover.getDouble(latitude, longitude));
    }

    double calcManningSurface(double volume) {
      if (volume < minToFlow)
        return 0;

      double height = volume / nodeArea;
      double vel = calcManning(.025, height, max(height / (2 * distAlong), slope));
      return min(vel, maxVelocity); //90 * sqrt(height / 2)); //(water terminal velocity)
    }
  };

  class HydroRiverNode : public HydroNode  {
  private:
    static const double minRiverToFlow = .001;
    static const double maxRiverVelocity = 4;

  public:
    HydroRiverNode(double nodeArea, double distAlong, double distAcross, double slope, Measure latitude, Measure longitude)
      : HydroNode(nodeArea, distAlong, distAcross, slope, minRiverToFlow, maxRiverVelocity, latitude, longitude) {
    }

    virtual HydroNode* clone(map<HydroNode*, HydroNode*>& translate) {
      HydroNode* copy = new HydroRiverNode(*this);
      copy->translateEdges(*this, translate);

      return copy;
    }

    virtual double calcVelocity() {
      return calcManningRiver(precipVolume + meltVolume);
    }

    virtual double calcEvaporationRate(GeographicMap<double>& temps, GeographicMap<double>& snowCover, GeographicMap<double>& relhums, GeographicMap<double>& lights, double fracyear, double surfaceEvaporationFactor, double riverEvaporationFactor) {
      return riverEvaporationFactor * HydroNode::calcEvaporation(temps.getDouble(latitude, longitude) - 273.15, relhums.getDouble(latitude, longitude), lights.getDouble(latitude, longitude), latitude, max(fracyear * 12 + 1, 12.0));
    }

    double calcRiverWidth(double volume) {
      return 2 * sqrt((3 / M_PI) * volume / distAlong);
    }

    double calcManningRiver(double volume) {
      if (volume < minToFlow)
        return 0;

      // 1/3 of circle calculation
      double width = calcRiverWidth(volume);
      if (width <= 0)
        return 0;

      double r, height;
      if (width > distAcross) {
        r = volume / nodeArea;
        height = r;
      } else {
        r = (volume / distAlong) / (M_PI * width / 3);
        height = width / 2;
      }

      double vel = calcManning(.033, r, max(height / (2 * distAlong), slope));
      return min(vel, maxVelocity); //90 * sqrt(height)); // terminal velocity
    }
  };

  class HydroOutputNode : public HydroNode  {
  public:
    HydroOutputNode()
      : HydroNode(0, 0, 0, 0, 0, 0, Measure(0, Inds::lat), Measure(0, Inds::lon)) {
    }

    virtual HydroNode* clone(map<HydroNode*, HydroNode*>& translate) {
      HydroNode* copy = new HydroOutputNode(*this);
      copy->translateEdges(*this, translate);

      return copy;
    }

    virtual void stepModel(double dt) {
      // ignore
    }

    virtual void stepPost() {
      // ignore
    }

    virtual double calcVelocity() {
      return 0;
    }

    virtual double calcEvaporationRate(GeographicMap<double>& temps, GeographicMap<double>& snowCover, GeographicMap<double>& relhums, GeographicMap<double>& lights, double fracyear, double surfaceEvaporationFactor, double riverEvaporationFactor) {
      return 0;
    }

    void reset() {
      precipVolume = precipVolumeAfter;
      meltVolume = meltVolumeAfter;
      if (confWeightAfter == 0)
        conf = 0;
      else
        conf = confAfter / confWeightAfter;

      precipVolumeAfter = meltVolumeAfter = confAfter = confWeightAfter = 0;
    }

    // Serialization

    static HydroOutputNode* streamExtractPointer(istream& in, PointerReference& reference) {
      throw logic_error("HydroOutputNode::streamExtractPointer not implemented yet!");
    }

  };

  class HydroNet {
  protected:
    GeographicMap<double> mask_coarse;
    map<unsigned, HydroNode*> surfaces_coarse; // mask_coarse index -> surface node
    map<unsigned, list<HydroNode*>*> nodemap_coarse; // mask_coarse index -> many nodes

  public:
    HydroNet(GeographicMap<double>& mask_coarse)
      : mask_coarse(mask_coarse) {
    }

    HydroNet(HydroNet& copy, map<HydroNode*, HydroNode*>& translate)
      : mask_coarse(copy.mask_coarse) {

      for (map<unsigned, HydroNode*>::iterator it = copy.surfaces_coarse.begin(); it != copy.surfaces_coarse.end(); it++) {
        surfaces_coarse[it->first] = HydroNode::getCopy(it->second, translate);
      }
      for (map<unsigned, list<HydroNode*>*>::iterator it = copy.nodemap_coarse.begin(); it != copy.nodemap_coarse.end(); it++) {
        list<HydroNode*>* nodes = new list<HydroNode*>();
        for (list<HydroNode*>::iterator lit = it->second->begin(); lit != it->second->end(); lit++)
          nodes->push_back(HydroNode::getCopy(*lit, translate));

        nodemap_coarse[it->first] = nodes;
      }
    }

    virtual ~HydroNet() {
      throw runtime_error("Deleting HydroNets is not currently implemented.");
    }

    HydroOutputNode* test(DInfinityMap& direction_fine) {
      return NULL;
    }

    HydroOutputNode* generate(DInfinityMap& direction_fine, GeographicMap<bool>& mask_fine, GeographicMap<double>& slope_fine, double min_dist) {
      HydroOutputNode* out = new HydroOutputNode();
      HydroOutputNode* ignore = new HydroOutputNode();

      map<unsigned, HydroNode*> rivers;

      surfaces_coarse.clear();
      nodemap_coarse.clear();

      // Loop over all surface cells
      for (unsigned rr = 0; rr < mask_coarse.getLatitudes().count(); rr++)
        for (unsigned cc = 0; cc < mask_coarse.getLongitudes().count(); cc++) {
          cout << "Generating cell " << rr << ", " << cc << endl;

          double total_dist = 0, total_slope = 0, weight_slope = 0;
          unsigned paths = 0;

          map<unsigned, pair<double, HydroNode*> > edges;

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
                  // connect this to the river out
                  if (rivers.find(mask_fine.getIndex(ll.first, ll.second)) == rivers.end())
                    generateRiver(ll.first, ll.second, rivers, direction_fine, mask_fine, slope_fine, out, ignore, min_dist);

                  addEdge(edges, mask_fine.getIndex(ll.first, ll.second), llp.second, rivers[mask_fine.getIndex(ll.first, ll.second)]);
                } else {
                  total_dist += llp.second * mask_fine.calcCellSpan(lat, lon);
                  total_slope += llp.second * slope_fine.getDouble(lat, lon);
                  weight_slope += llp.second;

                  followRiver(ll.first, ll.second, llp.second, pending, edges, direction_fine, mask_fine, out, ignore);
                }

                pending.pop();
              }
            }

          Measure latitude(Inds::lat), longitude(Inds::lon);
          mask_coarse.calcLatitudeLongitude(rr, cc, latitude, longitude);

          HydroNode* surface = new HydroSurfaceNode(mask_coarse.calcArea(rr, cc) * mask_coarse.getCellConst(rr, cc),
                                                    max(total_dist / paths, min_dist / 4), mask_coarse.calcDistance(rr - 1, cc - 1, rr + 1, cc + 1) / (2 * M_SQRT2),
                                                    total_slope / weight_slope, latitude, longitude);

          // divide edges weights by #paths
          surface->setEdges(edges, paths);
          surfaces_coarse[mask_coarse.getIndex(rr, cc)] = surface;
          addToNodeMap(surface, mask_coarse.getIndex(rr, cc));
        }

      return out;
    }

    // generate a river out of this cell
    void generateRiver(Measure lat, Measure lon, map<unsigned, HydroNode*>& rivers, DInfinityMap& direction_fine, GeographicMap<bool>& mask_fine, GeographicMap<double>& slope_fine, HydroOutputNode* out, HydroOutputNode* ignore, double min_dist, double bigger = 0) {
      //cout << "Generating river at " << lat << ", " << lon << " (" << bigger << ")" << endl;
      rivers[mask_fine.getIndex(lat, lon)] = ignore; // temporary

      map<unsigned, pair<double, HydroNode*> > edges;
      double total_dist = 0, total_slope = 0, weight_slope = 0;

      unsigned rr = mask_coarse.getLatitudes().inRange(lat), cc = mask_coarse.getLongitudes().inRange(lon);
      bool offedge = false;

      // Follow this out of the cell
      queue< pair< pair<Measure, Measure>, double> > pending;
      pending.push(pair< pair<Measure, Measure>, double>(pair<Measure, Measure>(lat, lon), 1.0));

      while (!pending.empty()) {
        pair< pair<Measure, Measure>, double> llp = pending.front();
        pair<Measure, Measure> ll = llp.first;
        if (ll.first < mask_coarse.getLatitudes().getCellMin(rr) - mask_coarse.getLatitudes().getWidths() * bigger ||
            ll.first >= mask_coarse.getLatitudes().getCellMax(rr) + mask_coarse.getLatitudes().getWidths() * bigger ||
            ll.second < mask_coarse.getLongitudes().getCellMin(cc) - mask_coarse.getLongitudes().getWidths() * bigger ||
            ll.second >= mask_coarse.getLongitudes().getCellMax(cc) + mask_coarse.getLongitudes().getWidths() * bigger) {
          //cout << "Exit at " << ll.first << ", " << ll.second << endl;
          // is this the output node?
          if (mask_fine.getDouble(ll.first, ll.second) == 0)
            rivers[mask_fine.getIndex(ll.first, ll.second)] = out;

          if (rivers.find(mask_fine.getIndex(ll.first, ll.second)) == rivers.end())
            generateRiver(ll.first, ll.second, rivers, direction_fine, mask_fine, slope_fine, out, ignore, min_dist);

          addEdge(edges, mask_fine.getIndex(ll.first, ll.second), llp.second, rivers[mask_fine.getIndex(ll.first, ll.second)]);
        } else {
          //cout << "Cont to " << ll.first << ", " << ll.second << endl;
          total_dist += llp.second * mask_fine.calcCellSpan(lat, lon);
          total_slope += llp.second * slope_fine.getDouble(lat, lon);
          weight_slope += llp.second;

          bool riverOffedge = followRiver(ll.first, ll.second, llp.second, pending, edges, direction_fine, mask_fine, out, ignore);
          offedge = offedge || riverOffedge;
        }

        pending.pop();
      }

      if (total_dist < min_dist && !offedge && bigger < 2) {
        //cout << "Total Dist: " << total_dist << endl;
        return generateRiver(lat, lon, rivers, direction_fine, mask_fine, slope_fine, out, ignore, min_dist, bigger + 1);
      }
      if (total_dist < min_dist)
        total_dist = min_dist;

      // actually make the river and add it to rivers
      HydroNode* river = new HydroRiverNode(total_dist * mask_fine.calcCellSpan(lat, lon),
                                            total_dist, mask_fine.calcCellSpan(lat, lon),
                                            total_slope / weight_slope, lat, lon);

      river->setEdges(edges);
      rivers[mask_fine.getIndex(lat, lon)] = river;
      addToNodeMap(river, mask_coarse.getIndex(lat, lon));
    }

    void addEdge(map<unsigned, pair<double, HydroNode*> >& edges, unsigned index, double weight, HydroNode* target) {
      if (edges.find(index) == edges.end())
        edges[index] = pair<double, HydroNode*>(weight, target);
      else
        edges[index].first += weight;
    }

    // returns true if flows off edge
    bool followRiver(Measure lat, Measure lon, double weight, queue< pair< pair<Measure, Measure>, double> >& pending, map<unsigned, pair<double, HydroNode*> >& edges, DInfinityMap& direction_fine, GeographicMap<bool>& mask_fine, HydroOutputNode* out, HydroOutputNode* ignore) {
      if (weight < .01) {
        addEdge(edges, mask_fine.getIndex(lat, lon), weight, ignore);
        return false;
      }

      Measure lat0(Inds::lat), lon0(Inds::lon), lat1(Inds::lat), lon1(Inds::lon);
      double portion0;

      bool offedge = false;

      direction_fine.getDirections(lat, lon, lat0, lon0, lat1, lon1, portion0);
      if (portion0 > 0) {
        if (mask_fine.getDouble(lat0, lon0) == 0) {
          if (portion0 < 1)
            addEdge(edges, mask_fine.getIndex(lat, lon), portion0 * weight, ignore);
          else
            addEdge(edges, mask_fine.getIndex(lat, lon), portion0 * weight, out);
          offedge = true;
        } else {
          if (lat0 >= mask_fine.getLatitudes().getMin() + mask_fine.getLatitudes().getWidths() &&
              lat0 < mask_fine.getLatitudes().getMax() - mask_fine.getLatitudes().getWidths() &&
              lon0 >= mask_fine.getLongitudes().getMin() + mask_fine.getLongitudes().getWidths() &&
              lon0 < mask_fine.getLongitudes().getMax() - mask_fine.getLongitudes().getWidths())
            pending.push(pair< pair<Measure, Measure>, double>(pair<Measure, Measure>(lat0, lon0), portion0 * weight));
          else {
            addEdge(edges, mask_fine.getIndex(lat, lon), portion0 * weight, out);
            offedge = true;
          }
        }
      }
      if (portion0 < 1) {
        if (mask_fine.getDouble(lat1, lon1) == 0) {
          if (portion0 > 0)
            addEdge(edges, mask_fine.getIndex(lat, lon), (1 - portion0) * weight, ignore);
          else
            addEdge(edges, mask_fine.getIndex(lat, lon), (1 - portion0) * weight, out);
          offedge = true;
        } else {
          if (lat1 >= mask_fine.getLatitudes().getMin() + mask_fine.getLatitudes().getWidths() &&
              lat1 < mask_fine.getLatitudes().getMax() - mask_fine.getLatitudes().getWidths() &&
              lon1 >= mask_fine.getLongitudes().getMin() + mask_fine.getLongitudes().getWidths() &&
              lon1 < mask_fine.getLongitudes().getMax() - mask_fine.getLongitudes().getWidths())
            pending.push(pair< pair<Measure, Measure>, double>(pair<Measure, Measure>(lat1, lon1), (1 - portion0) * weight));
          else {
            addEdge(edges, mask_fine.getIndex(lat, lon), (1 - portion0) * weight, out);
            offedge = true;
          }
        }
      }

      return offedge;
    }

    void addToNodeMap(HydroNode* node, unsigned index) {
      list<HydroNode*>* nodes = NULL;

      if (nodemap_coarse.find(index) == nodemap_coarse.end()) {
        nodes = new list<HydroNode*>();
        nodemap_coarse[index] = nodes;
      } else
        nodes = nodemap_coarse[index];

      nodes->push_back(node);
    }

    void sumNodeMapVolumes(unsigned rr, unsigned cc, double& precipVolume, double& meltVolume) {
      precipVolume = meltVolume = 0;
      list<HydroNode*>* nodes = NULL;
      unsigned index = mask_coarse.getIndex(rr, cc);

      if (nodemap_coarse.find(index) == nodemap_coarse.end())
        return;

      nodes = nodemap_coarse[index];
      for (list<HydroNode*>::iterator it = nodes->begin(); it != nodes->end(); it++) {
        precipVolume += (*it)->getPrecipVolume();
        meltVolume += (*it)->getMeltVolume();
      }
    }

    double calculateMaximumStep() {
      double maxstep = FLT_MAX;

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

    void updateDay(GeographicMap<double>& temps, GeographicMap<double>& snowCover, GeographicMap<double>& precipChanges, double fracyear, double surfaceEvaporationFactor, double riverEvaporationFactor) {
      // Construct relhums and lights
      GeographicMap<double>& relhums = .3 * (precipChanges > 0) + .5; // raining air is 60-100% relhum.
      GeographicMap<double>& lights = (-0.5 * (precipChanges > 0) + 0.5);

      map<unsigned, HydroNode*>::iterator it;
      for (it = surfaces_coarse.begin(); it != surfaces_coarse.end(); it++)
        it->second->resetRecursiveCalculation();

      for (it = surfaces_coarse.begin(); it != surfaces_coarse.end(); it++)
        it->second->recursiveUpdateEvaporationRate(temps, snowCover, relhums, lights, fracyear, surfaceEvaporationFactor, riverEvaporationFactor);
    }

    void step(double dt, GeographicMap<double>& precipChanges, GeographicMap<double>& meltChanges, GeographicMap<double>& changesConf) {
      map<unsigned, HydroNode*>::iterator it;
      for (it = surfaces_coarse.begin(); it != surfaces_coarse.end(); it++)
        it->second->stepModel(dt);

      for (unsigned rr = 0; rr < mask_coarse.getLatitudes().count(); rr++)
        for (unsigned cc = 0; cc < mask_coarse.getLongitudes().count(); cc++) {
          double fraction = mask_coarse.getCellConst(rr, cc);
          if (fraction == 0)
            continue;

          double precip = precipChanges.getCellConst(rr, cc) * fraction * dt / DividedRange::toTimespan(1).getValue();
          double melt = meltChanges.getCellConst(rr, cc) * fraction * dt / DividedRange::toTimespan(1).getValue();
          double conf = changesConf.getCellConst(rr, cc);
          surfaces_coarse[mask_coarse.getIndex(rr, cc)]->addStepVolume(precip, melt, conf, 1.0);
        }

      for (it = surfaces_coarse.begin(); it != surfaces_coarse.end(); it++)
        it->second->stepPost();
    }

    // diagnostics

    list<pair<pair<HydroNode*, HydroNode*>, double> > getAllEdges() {
      list<pair<pair<HydroNode*, HydroNode*>, double> > edges;

      map<unsigned, HydroNode*>::iterator it;
      for (it = surfaces_coarse.begin(); it != surfaces_coarse.end(); it++)
        it->second->resetRecursiveCalculation();

      for (it = surfaces_coarse.begin(); it != surfaces_coarse.end(); it++) {
        list<pair<pair<HydroNode*, HydroNode*>, double> > edgeSet = it->second->recursiveGetAllEdges();
        edges.insert(edges.end(), edgeSet.begin(), edgeSet.end());
      }

      return edges;
    }

    // Serializable protocol

    friend void* nodeListConstructor(istream& in, PointerReference& reference);

    virtual istream& streamExtract(istream& in, PointerReference& reference) {
      unsigned surfaces_count;
      in >> surfaces_count;

      for (unsigned ii = 0; ii < surfaces_count; ii++) {
	unsigned index;
	in >> index;
	throw runtime_error("Need to have a extractor dictionary, and save the class names!");
	//XXX: HydroNode* node = HydroNode::streamExtractPointer(in, reference);
	//surfaces_coarse.insert(pair<unsigned, HydroNode*>(index, node));
      }

      unsigned nodemap_count;
      in >> nodemap_count;

      for (unsigned ii = 0; ii < nodemap_count; ii++) {
	unsigned index;
	in >> index;

	list<HydroNode*>* nodeList = (list<openworld::HydroNode*>*) reference.streamExtractPointer(in, nodeListConstructor);
	nodemap_coarse.insert(pair<unsigned, list<HydroNode*>*>(index, nodeList));
      }

      return in;
    }

    virtual ostream& streamInsert(ostream& os, PointerTracker& tracker) const {
      os << surfaces_coarse.size() << " ";

      for (map<unsigned, HydroNode*>::const_iterator it = surfaces_coarse.begin() ; it != surfaces_coarse.end(); it++) {
	os << it->first << " ";
	it->second->streamInsertPointer(os, tracker);
      }

      os << nodemap_coarse.size() << " ";

      for (map<unsigned, list<HydroNode*>*>::const_iterator it = nodemap_coarse.begin() ; it != nodemap_coarse.end(); it++) {
	os << it->first << " ";
	if (tracker.streamInsertPointer(os, it->second)) {
	  os << it->second->size() << " ";

	  for (list<HydroNode*>::iterator jt = it->second->begin(); jt != it->second->end(); ++jt)
	    (*jt)->streamInsertPointer(os, tracker);
	}
      }

      return os;
    }

  };
}

#endif
