#ifndef EVENT_DISPLAY_H
#define EVENT_DISPLAY_H

#include <array>
#include <vector>

class TApplication;

namespace cluster_utils {
struct LayerClusters;
}

void runPerEventDisplay(
    const std::vector<std::array<cluster_utils::LayerClusters, 3>>& clustersPerEvent,
    bool showDisplay,
    TApplication* appPtr,
    double layerSpacingCm);

#endif  // EVENT_DISPLAY_H
