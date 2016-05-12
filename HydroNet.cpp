#include "HydroNet.h"

namespace openworld {
  const double HydroSurfaceNode::minSurfaceToFlow = .1;
  const double HydroSurfaceNode::maxSurfaceVelocity = 1;
  const double HydroRiverNode::minRiverToFlow = .001;
  const double HydroRiverNode::maxRiverVelocity = 4;

  void* nodeListConstructor(istream& in, PointerReference& reference) {
    unsigned count;
    in >> count;

    list<HydroNode*>* nodes = new list<HydroNode*>();
    for (unsigned ii = 0; ii < count; ii++) {
      HydroNode* node = HydroNode::streamExtractPointer(in, reference);
      nodes->push_back(node);
    }

    return nodes;
  }
}
