#include "HydroNet.h"

namespace openworld {
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
